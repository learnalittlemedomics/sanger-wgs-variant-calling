#!/usr/bin/env python3
"""Mock integration: branching, official split indices, strict join, empty output."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

DRIVER = Path(__file__).resolve().parents[2]/'scripts/sanger-split/brass-scatter.py'
spec=importlib.util.spec_from_file_location('brass',str(DRIVER))
brass=importlib.util.module_from_spec(spec)
spec.loader.exec_module(brass)

MOCK = r'''#!/usr/bin/env python3
import os,sys
from pathlib import Path
args=sys.argv[1:]
stage=args[args.index('-p')+1]
out=Path(args[args.index('-o')+1]); tmp=out/'tmpBrass'
tmp.mkdir(parents=True,exist_ok=True)
def put(p,s='x'):
 p.parent.mkdir(parents=True,exist_ok=True); p.write_text(s)
if stage in ('normcn','filter','grass','tabix'):
 assert Path(args[args.index('-ss')+1]).is_file(), 'Missing ASCAT'
if stage=='input':
 assert 'brass_input_lib' in os.environ['PERL5LIB']
 role=int(args[args.index('-i')+1]); assert args[args.index('-c')+1]=='30'
 put(tmp/('T.brm.bam' if role==1 else 'N.brm.bam'))
elif stage=='cover':
 assert not (tmp/'T.brm.bam').exists(), 'Cover was incorrectly gated by input'
 put(tmp/'cover/raw')
elif stage=='group':
 assert (tmp/'T.brm.bam').is_file(); put(out/'groups')
elif stage=='isize': put(out/'insert_size')
elif stage=='merge':
 assert (tmp/'cover/raw').is_file(); put(tmp/'cover/merged')
elif stage=='normcn':
 assert (tmp/'cover/merged').is_file(); put(out/'cn')
elif stage=='filter':
 for x in ('cn','groups','insert_size'): assert (out/x).is_file(), x
 put(out/'groups.clean.bedpe')
elif stage=='split':
 for i in range(int(os.environ['TEST_SPLITS'])): put(tmp/'split'/('split.%07d'%i),'')
elif stage=='assemble':
 i=int(args[args.index('-i')+1]); assert args[args.index('-c')+1]=='1'
 assert (tmp/'split'/('split.%07d'%(i-1))).is_file()
 put(tmp/'assemble'/('bedpe.%07d'%(i-1)), '' if os.environ['TEST_SPLITS']=='1' else str(i))
elif stage=='grass':
 assert len(list((tmp/'assemble').glob('bedpe.*')))==int(os.environ['TEST_SPLITS'])
 put(out/'final.vcf')
elif stage=='tabix': assert (out/'final.vcf').is_file()
else: raise RuntimeError(stage)
put(tmp/'progress'/(stage+('.'+args[args.index('-i')+1] if '-i' in args else '')))
'''

class Integration(unittest.TestCase):
 def test_chain_and_empty(self):
  for count in (7,1):
   with tempfile.TemporaryDirectory() as d:
    root=Path(d); mock=root/'mock.py'; mock.write_text(MOCK); mock.chmod(0o755)
    lib=root/'lib/Sanger/CGP/Brass'; lib.mkdir(parents=True)
    (lib/'Implement.pm').write_text('package Sanger::CGP::Brass::Implement;\nmy $helper_threads=15;\n$helper_threads = 6 if($helper_threads > 6);\n1;\n')
    env=os.environ.copy(); env.update(PERL5LIB=str(root/'lib'),TEST_SPLITS=str(count))
    archives={}
    original=os.getcwd()
    try:
     ascat=root/'ascat'; ascat.mkdir(); os.chdir(str(ascat))
     Path('ascat.txt').write_text('rho 0.7')
     brass.pack('ascat.tar',['ascat.txt']); asc=ascat/'ascat.tar'
     def task(stage,deps=(),group=0):
      work=root/(stage+str(group)); work.mkdir(); os.chdir(str(work))
      for arc in deps: brass.extract(arc)
      cmd=[sys.executable,str(DRIVER.with_name('brass-gather.py') if stage=='gather' else DRIVER),'run','--stage',stage,'--cpus','30','--command',
           str(mock)+' -j 4 -k 4 -c 30 -o WGS_T_vs_N/brass -ss ascat.txt',
           '--brass-dir','WGS_T_vs_N/brass','--group',str(group),'--workers','6']
      subprocess.check_call(cmd,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
      return work
     a=task('input',group=1); b=task('input',group=2); c=task('cover'); i=task('isize')
     g=task('group',[a/'BRASS_input_1.stage.tar',b/'BRASS_input_2.stage.tar'])
     n=task('normcn',[c/'BRASS_cover.stage.tar',asc])
     p=task('plan',[g/'BRASS_group.stage.tar',n/'BRASS_normcn.stage.tar',i/'BRASS_isize.stage.tar',asc])
     plan=json.loads((p/'brass-plan.json').read_text())
     parts=[]
     for index in range(1,len(plan['groups'])+1):
      w=task('assemble',[p/'brass-seed.tar'],index)
      parts.append(w/('brass-part-%d.tar'%index))
     final=task('gather',parts+[p/'BRASS_plan.stage.tar',asc])
     self.assertTrue((final/'BRASS.stage.tar.gz').is_file())
     if len(parts)>1:
      with self.assertRaises(subprocess.CalledProcessError):
       # Missing group must fail before annotation.
       task('gather',parts[:1]+[p/'BRASS_plan.stage.tar',asc],99)
     os.chdir(str(root)); Path('conflict').write_text('a'); brass.pack('a.tar',['conflict'])
     Path('conflict').write_text('b')
     with self.assertRaises(RuntimeError): brass.extract('a.tar')
    finally: os.chdir(original)

if __name__=='__main__': unittest.main()
