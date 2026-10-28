#!/usr/bin/env python3
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]
WORKER=ROOT/'scripts/sanger-split/pindel-scatter.py'
MOCK=r'''#!/usr/bin/env python3
import gzip,json,os,sys
from pathlib import Path
args=sys.argv[1:]
def arg(k):return args[args.index(k)+1]
stage=arg('-p');out=Path(arg('-o'));tmp=out/'tmpPindel'
idx=int(arg('-i')) if '-i' in args else 0
if os.environ.get('FAIL_PINDEL')=='oom' and stage=='pindel':
 print('cannot allocate memory');sys.exit(1)
prefix='Sanger_CGP_Pindel_Implement_'
if stage=='input':
 assert arg('-c')==os.environ.get('EXPECTED_INPUT_THREADS','4')
 sd=tmp/('T' if idx==1 else 'N');sd.mkdir(exist_ok=True)
 for chrom in ('chr1','chr2'):
  with gzip.open(sd/(chrom+'.txt.gz'),'wt') as f:f.write('reads')
elif stage=='pindel':
 assert arg('-c')=='1'
 for sample in ('T','N'):
  assert sorted(p.name for p in (tmp/sample).glob('*.txt.gz'))==['chr1.txt.gz','chr2.txt.gz']
 chrom=['chr1','chr2'][idx-1];d=tmp/'pout';d.mkdir(exist_ok=True)
 for ext in ('D','SI'):(d/(chrom+'_'+chrom+'_'+ext)).write_text('calls')
elif stage=='pin2vcf':
 assert (tmp/'progress'/(prefix+'pindel.'+str(idx))).exists()
 chrom=['chr1','chr2'][idx-1];d=tmp/'vcf';d.mkdir(exist_ok=True)
 (d/(chrom+'_pindel.vcf')).write_text('##fileformat=VCFv4.2\n')
 (d/(chrom+'_pindel.sam')).write_text('alignment')
elif stage=='merge':
 assert len(list((tmp/'vcf').glob('*.vcf')))==2
 (out/'T_vs_N.vcf.gz').write_bytes(b'merged')
elif stage=='flag':
 assert (out/'T_vs_N.vcf.gz').exists()
 with gzip.open(out/'T_vs_N.flagged.vcf.gz','wt') as f:f.write('##fileformat=VCFv4.2\n')
 (out/'T_vs_N.flagged.vcf.gz.tbi').write_text('index')
 (out/'T_vs_N.germline.bed').write_text('chr1\t0\t10\n')
else:raise ValueError(stage)
method={'pin2vcf':'pindel_to_vcf','merge':'merge_and_bam'}.get(stage,stage)
(tmp/'progress'/(prefix+method+'.'+str(idx))).touch()
'''
PERL=r'''#!/usr/bin/env python3
import json,sys
from pathlib import Path
o=json.load(sys.stdin)
r=dict(version='mock',tumour_name='T',normal_name='N',tumour_dir='T',normal_dir='N')
if o['layout']:
 seq=sorted(set(p.name[:-7] for s in ('T','N') for p in (Path(o['tmp'])/s).glob('*.txt.gz')))
 r.update(jobs=len(seq),chromosomes=seq)
print(json.dumps(r))
'''
class TestPindel(unittest.TestCase):
 def test_official_stages_and_rejection(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);bins=root/'bin';bins.mkdir()
   for name,body in [('pindel.pl',MOCK),('perl',PERL),('bgzip', '#!/usr/bin/env python3\nimport gzip,sys\nsys.stdout.buffer.write(gzip.compress(open(sys.argv[-1],"rb").read()))\n'),('tabix','#!/usr/bin/env python3\nimport sys\nopen(sys.argv[-1]+".tbi","w").write("index")\n')]:
    p=bins/name;p.write_text(body);p.chmod(0o755)
   env=dict(os.environ,PATH=str(bins)+os.pathsep+os.environ['PATH'])
   parts=[]
   def run(mode,name,extra=(),code=0,oom=False,worker=WORKER):
    work=root/name;work.mkdir(exist_ok=True)
    (work/'split_scripts').symlink_to(ROOT/'scripts/sanger-split') if not (work/'split_scripts').exists() else None
    cmd='pindel.pl -o '+str(work/'WGS_T_vs_N/pindel')+' -t t.cram -n n.cram -r genome.fa -e chrM -c 1'
    r=subprocess.run([sys.executable,str(worker),mode,'--command',cmd]+list(extra),cwd=str(work),env=dict(env,FAIL_PINDEL='oom' if oom else '',EXPECTED_INPUT_THREADS='16' if worker.name=='pindel-input16.py' else '4'),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    self.assertEqual(r.returncode,code,r.stdout+'\n'+r.stderr)
    return work
   for i in (1,2):
    work=run('input','input'+str(i),['--index',str(i)])
    parts.append(work/('pindel_input_'+str(i)))
   for i in (1,2):
    newer=run('input','input16_'+str(i),['--index',str(i)],worker=ROOT/'scripts/sanger-split/pindel-input16.py')
    self.assertEqual(json.loads((newer/('pindel_input_'+str(i))/'report.json').read_text()),json.loads((parts[i-1]/'report.json').read_text()))
   common=['--inputs']+[str(p) for p in parts]
   plan=run('plan','plan',common)/'pindel-plan.json'
   self.assertEqual(json.loads(plan.read_text())['chromosomes'],['chr1','chr2'])
   outputs=[]
   for i in (1,2):
    args=['--index',str(i),'--plan',str(plan)]+common
    if i==1:run('chromosome','chr'+str(i),args,137,True)
    work=run('chromosome','chr'+str(i),args)
    outputs.append(work/('pindel_part_'+str(i)))
   run('gather','missing',['--plan',str(plan),'--parts',str(outputs[0])],1)
   run('gather','duplicate',['--plan',str(plan),'--parts',str(outputs[0]),str(outputs[0])],1)
   work=run('gather','gather',['--plan',str(plan),'--parts']+[str(p) for p in outputs])
   with tarfile.open(str(work/'cgpPindel.stage.tar.gz')) as tf:
    self.assertIn('WGS_T_vs_N/pindel/T_vs_N.flagged.vcf.gz',tf.getnames())
    self.assertIn('WGS_T_vs_N/pindel/T_vs_N.germline.bed.gz.tbi',tf.getnames())
   final=run('gather','final',['--plan',str(plan),'--parts']+[str(p) for p in outputs],worker=ROOT/'scripts/sanger-split/pindel-finalize.py')
   self.assertFalse((final/'cgpPindel.stage.tar.gz').exists())
   with tarfile.open(str(final/'pindel-germline.stage.tar.gz')) as tf:
    self.assertEqual(sorted(tf.getnames()),['WGS_T_vs_N/pindel/T_vs_N.germline.bed.gz','WGS_T_vs_N/pindel/T_vs_N.germline.bed.gz.tbi'])
   for name in ('T_vs_N.flagged.vcf.gz','T_vs_N.germline.bed.gz'):
    self.assertEqual(gzip.open(str(work/'WGS_T_vs_N/pindel'/name),'rb').read(),gzip.open(str(final/'WGS_T_vs_N/pindel'/name),'rb').read())
   pack=root/'pack';pack.mkdir()
   for name in ('WGS_T_vs_N','timings'):(pack/name).symlink_to(final/name,target_is_directory=True)
   nf=(ROOT/'modules/sanger-split.nf').read_text()
   script=nf.split('process PindelPack {',1)[1].split('script:',1)[1].split('"""',2)[1]
   script=script.replace('${pair}','WGS_T_vs_N').replace('${timings}','timings').replace('\\$','$')
   subprocess.run(['bash','-euc',script],cwd=str(pack),check=True)
   with tarfile.open(str(work/'cgpPindel.stage.tar.gz')) as old, tarfile.open(str(pack/'cgpPindel.stage.tar.gz')) as new:
    # Parent directory entry differs, result files must not.
    self.assertEqual(sorted(m.name for m in old if m.isfile()),sorted(m.name for m in new if m.isfile()))
   self.assertTrue((parts[0]/'data/chr1.txt.gz').is_file())
if __name__=='__main__':unittest.main(verbosity=2)
