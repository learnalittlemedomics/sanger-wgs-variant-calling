#!/usr/bin/env python3
"""Integration tests with mock CaVEMan; never submits a Slurm job or reads real CRAM.
Run: python3 tests/scatter/test_scatter.py
"""
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
WORKER = ROOT/'scripts/sanger-split/caveman-scatter.py'
spec = importlib.util.spec_from_file_location('scatter', WORKER)
scatter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scatter)

MOCK = r'''#!/usr/bin/env python3
import gzip,json,os,sys,time,fcntl,atexit
from pathlib import Path
args=sys.argv[1:]
def arg(k): return args[args.index(k)+1]
stage=arg('-p')
if stage in ('mstep','estep'):
 tracker=os.environ['MOCK_TRACKER']
 def track(delta):
  with open(tracker,'a+') as f:
   fcntl.flock(f,fcntl.LOCK_EX);f.seek(0)
   state=json.loads(f.read() or '{"active":0,"peak":0}')
   state['active']+=delta;state['peak']=max(state['peak'],state['active'])
   f.seek(0);f.truncate();json.dump(state,f);f.flush()
 # Deliberately stagger startup beyond the old 0.1-second overlap window.
 time.sleep(float(os.environ.get('MOCK_STARTUP_STAGGER', '0')) * int(arg('-i')))
 track(1);atexit.register(track,-1)
 # First wave rendezvous: test configured concurrency without depending on
 # process launch speed or filesystem latency. Later waves pass immediately.
 barrier=Path(os.environ['MOCK_BARRIER'])
 (barrier/(str(os.getpid())+'.ready')).touch()
 target=int(os.environ['MOCK_BARRIER_TARGET'])
 deadline=time.monotonic()+60
 while len(list(barrier.glob('*.ready'))) < target:
  if time.monotonic() >= deadline:
   print('Mock concurrency rendezvous timed out: expected '+str(target),flush=True)
   sys.exit(1)
  time.sleep(0.02)
if os.environ.get('MOCK_FAIL')=='oom':
 print('ERROR: cannot allocate memory',flush=True);sys.exit(1)
if os.environ.get('MOCK_FAIL')=='ordinary':
 print('ERROR: invalid input',flush=True);sys.exit(1)
out=Path(arg('-o'));tmp=out/'tmpCaveman'
if stage=='split_concat':
 assert arg('-i')=='1'
 (tmp/'splitList').write_text((tmp/'split-fixture').read_text())
 (tmp/'progress'/'Sanger_CGP_Caveman_Implement_caveman_split_concat').touch()
 sys.exit(0)
rows=[line.split() for line in (tmp/'splitList').read_text().splitlines()]
assert int(arg('-t'))>=1
for key in ('-tb','-nb','-tc','-nc'):
 assert Path(arg(key)).is_file(),(key,arg(key))
prefix='Sanger_CGP_Caveman_Implement_caveman_'
(tmp/'progress').mkdir(exist_ok=True)
if stage in ('mstep','estep'):
 assert arg('-t')=='1'
 i=int(arg('-i'));g=int(arg('-l'))
 for n,(c,s,e) in enumerate(rows,1):
  if (n-1)%g!=i-1: continue
  p=tmp/'results'/c;p.mkdir(parents=True,exist_ok=True)
  stem=str(int(s)+1)+'_'+e
  if stage=='mstep':
   with gzip.open(p/(stem+'.covs'),'wb') as f:f.write(str(n).encode())
  else:
   assert (tmp/'cov_arr').read_bytes()==b'new model'
   for ext in ('muts.vcf.gz','snps.vcf.gz'):
    with gzip.open(p/(stem+'.'+ext),'wt') as f:f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+c+'\t'+str(int(s)+1)+'\t.\tA\tT\t.\t.\t.\n')
   (p/(stem+'.no_analysis.bed')).write_text('')
  (tmp/'progress'/(prefix+stage+'.'+str(n))).touch()
elif stage=='merge':
 assert len(list((tmp/'results').rglob('*.covs')))==len(rows)
 (tmp/'cov_arr').write_bytes(b'new model');(tmp/'prob_arr').write_bytes(b'new probabilities')
elif stage=='merge_results':
 assert len(list((tmp/'results').rglob('*.muts.vcf.gz')))==len(rows)
 for ext in ('muts.vcf.gz','snps.vcf.gz'):
  with gzip.open(out/('T_vs_N.'+ext),'wt') as f:f.write('merged')
elif stage=='add_ids':
 assert (out/'T_vs_N.muts.vcf.gz').is_file()
 (out/'ids_done').touch()
elif stage=='flag':
 assert (out/'ids_done').is_file()
 assert Path(arg('-in')).is_file()
 assert len(list((tmp/'progress').glob(prefix+'estep.*')))==len(rows)
 shutil=None
 with gzip.open(out/'T_vs_N.flagged.muts.vcf.gz','wt') as f:f.write('flagged')
else: raise RuntimeError(stage)
(tmp/'progress'/(prefix+stage)).touch()
'''


class ScatterTest(unittest.TestCase):
    def test_python35_path_boundaries(self):
        import ast
        tree = ast.parse(WORKER.read_text())
        counts = {'copy2': 2, 'copytree': 2, 'rmtree': 1, 'move': 2}
        checked = 0
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
                continue
            owner = node.func.value
            if not isinstance(owner, ast.Name):
                continue
            count = counts.get(node.func.attr, 0) if owner.id == 'shutil' else 0
            if owner.id == 'gzip' and node.func.attr == 'open':
                count = 1
            for arg in node.args[:count]:
                self.assertIsInstance(arg, ast.Call)
                self.assertIsInstance(arg.func, ast.Name)
                self.assertEqual(arg.func.id, 'str')
                checked += 1
        self.assertGreaterEqual(checked, 10)

    def test_pindel_parallelism(self):
        command = ROOT/'scripts/sanger-split/command-cgpPindel.sh'
        env = dict(os.environ)
        for key in ('OUTPUT_DIR','PROTOCOL','NAME_MT','NAME_WT','REF_BASE',
                    'BAM_MT_TMP','BAM_WT_TMP','ASSEMBLY','SPECIES','CONTIG_EXCLUDE'):
            env[key] = 'fixture'
        for cpus, workers, expected in ((30,16,16),(30,12,12),(30,8,8),(4,16,4)):
            env.update(CPU=str(cpus), PINDEL_WORKERS=str(workers))
            result = subprocess.run(['bash','-eu','-c',
                'declare -A do_parallel; source "$1"; printf "%s" "${do_parallel[cgpPindel]}"',
                'test',str(command)], env=env, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, universal_newlines=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(' -c '+str(expected)+' ',result.stdout)

    def test_groups(self):
        self.assertEqual(scatter.group_count(4802, 3.328, 1.5), 178)
        self.assertEqual(scatter.group_count(4802, 15.83, 1.5), 845)
        self.assertEqual(-(-scatter.group_count(4802, 3.328, 1.5)//6),30)
        self.assertEqual(-(-scatter.group_count(4802, 15.83, 1.5)//6),141)
        self.assertEqual(scatter.group_count(1, 100, 1.5), 1)
        with self.assertRaises(ValueError): scatter.group_count(1, 0, 1.5)

    def test_production_chain_and_failures(self):
        for packing, workers in ((1,1), (2,2), (6,6), (6,2)):
            with self.subTest(packing=packing, workers=workers):
                self.run_chain(packing,workers)

    def run_chain(self, packing,workers):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td); work=base/'work'; work.mkdir(); store=base/'store';store.mkdir()
            bins=base/'bin';bins.mkdir()
            (bins/'caveman.pl').write_text(MOCK);(bins/'caveman.pl').chmod(0o755)
            (bins/'perl').write_text('#!/bin/bash\nif [[ "$*" == *"-MSanger::CGP::Caveman"* ]]; then printf 1.13.14; else exec '+shutil.which('perl')+' "$@"; fi\n')
            (bins/'perl').chmod(0o755)
            (bins/'samtools').write_text('#!/bin/bash\nif [[ "$*" == *t.cram* ]]; then printf "@RG\\tID:1\\tSM:T\\n"; else printf "@RG\\tID:1\\tSM:N\\n"; fi\n')
            (bins/'samtools').chmod(0o755)
            (bins/'time').write_text('#!/bin/bash\n[[ "$1" != -v ]] || shift\nexec "$@"\n')
            (bins/'time').chmod(0o755)
            mock_scripts=base/'scripts'
            shutil.copytree(ROOT/'scripts/sanger-split',mock_scripts)
            runner=mock_scripts/'run-stage.sh'
            runner.write_text(runner.read_text().replace('/usr/bin/time',str(bins/'time')))
            env=dict(os.environ,PATH=str(bins)+os.pathsep+os.environ['PATH'],MOCK_TRACKER=str(base/'tracker.json'))
            pair='WGS_T_vs_N';cave=work/pair/'caveman'
            def inputs():
                for f in ('t.cram','t.cram.crai','t.cram.bas','n.cram','n.cram.crai','n.cram.bas'):
                    (work/f).touch()
                (work/'reference_files').mkdir();(work/'reference_files/genome.fa').write_text('>chr1\nA\n')
                (work/'reference_files/genome.fa.fai').write_text('chr1\t1000\n')
                (work/'ref_cache').mkdir()
                (work/'split_scripts').symlink_to(mock_scripts,target_is_directory=True)
                (work/'run.params').write_text("\n".join([
                    "OUTPUT_DIR='{}'".format(work),"REF_BASE='{}'".format(work/'reference_files'),
                    "BAM_MT='t.cram'","IDX_MT='t.cram.crai'","BAM_WT='n.cram'","IDX_WT='n.cram.crai'",
                    "PINDEL_MAXCPU=6","SPECIES=human","ASSEMBLY=GRCh38","CAVESPLIT=800000",
                    "SNVFLAG=flag.ini","CONTIG_EXCLUDE=chrM","PROTOCOL=WGS"]))
            def run(args,code=0,extra=None):
                p=subprocess.run(args,cwd=work,env=dict(env,**(extra or {})),capture_output=True,text=True)
                self.assertEqual(p.returncode,code,p.stdout+'\n'+p.stderr)
                return p
            inputs()
            (cave/'tmpCaveman/progress').mkdir(parents=True)
            (cave/'tmpCaveman/alg_bean').write_text('algorithm')
            (cave/'tmpCaveman/caveman.cfg.ini').write_text('CWD={}\n'.format(work))
            (cave/'tmpCaveman/split-fixture').write_text(''.join('chr1 {} {}\n'.format(i*100,(i+1)*100) for i in range(12)))
            for stage in ('setup','split'):(cave/'tmpCaveman/progress'/(scatter.PREFIX+stage)).touch()
            ad=work/pair/'ascat';ad.mkdir()
            (ad/'T.copynumber.caveman.csv').write_text('1,chr1,1,1000,2,1,2\n')
            (ad/'T.samplestatistics.txt').write_text('rho\t0.7\n')
            # End-to-end preparation uses the real shell and the unchanged scientific commands.
            run(['bash','split_scripts/run-stage.sh','CaVEMan_scatter_prepare','1'],extra={
                'CAVE_WORKERS_PER_JOB':str(packing),'CAVE_TARGET_HOURS':'0.1','CAVE_MSTEP_MINUTES':'2','CAVE_ESTEP_MINUTES':'3'})
            shutil.copytree(work/'scatter_context',store/'context')
            plan=json.loads((store/'context/plan.json').read_text())
            self.assertEqual((plan['mstep_worker_groups'],plan['estep_worker_groups']),(4,6))
            self.assertEqual((plan['mstep_groups'],plan['estep_groups']),(-(-4//packing),-(-6//packing)))
            self.assertEqual(plan['command'][plan['command'].index('-k')+1],'0.3')
            for name in ('ascat',): shutil.copytree(work/pair/name,store/name)
            def reset():
                shutil.rmtree(work);work.mkdir();inputs()
                (work/'scatter_context').symlink_to(store/'context',target_is_directory=True)
            def compute(mode,index=1,parts=(),model=False,extra=None,code=0):
                args=['python3',str(WORKER),'run','--mode',mode,'--index',str(index),'--workers',str(workers)]
                if parts: args+=['--parts']+[str(x) for x in parts]
                if model: args+=['--model',str(store/'model')]
                test_env=dict(extra or {})
                if mode in ('mstep','estep'):
                    lanes=min(packing,plan[mode+'_worker_groups']-(index-1)*packing)
                    test_env.update(MOCK_BARRIER=tempfile.mkdtemp(prefix='barrier-',dir=base),
                                    MOCK_BARRIER_TARGET=str(min(workers,lanes)),
                                    MOCK_STARTUP_STAGGER='0.15' if packing==6 else '0')
                return run(args,code,test_env)
            reset()
            # OOM propagated even if the Perl wrapper exits 1; ordinary errors stay 1.
            compute('mstep',extra={'MOCK_FAIL':'oom'},code=137)
            compute('mstep',extra={'MOCK_FAIL':'ordinary'},code=1)
            # Reusing the same task directory after a failed attempt is supported.
            mparts=[]
            for i in range(1,plan['mstep_groups']+1):
                if i>1: reset()
                compute('mstep',i);part=store/('m'+str(i))
                shutil.copytree(work/'result',part);mparts.append(part)
            reset();compute('merge',parts=mparts)
            shutil.copytree(work/'result',store/'model')
            reset();compute('merge',parts=[store/'m1',store/'m1'],code=1)
            reset();compute('merge',parts=[],code=1)
            for i in range(1,plan['estep_groups']+1):
                reset();compute('estep',i,model=True);shutil.copytree(work/'result',store/('e'+str(i)))
            reset();compute('gather',parts=[store/('e'+str(i)) for i in range(1,plan['estep_groups']+1)],model=True)
            archive=store/'CaVEMan.stage.tar.gz';shutil.copy2(work/archive.name,archive)
            with tarfile.open(archive) as tar:
                self.assertIn(pair+'/caveman/ids_done',tar.getnames())
            # Test the handoff into the ORIGINAL flag stage, not a new flag implementation.
            reset()
            subprocess.run(['tar','-xzf',str(archive),'-C',str(work)],check=True)
            shutil.copytree(store/'ascat',work/pair/'ascat')
            pindel=work/pair/'pindel';pindel.mkdir()
            (pindel/'T_vs_N.germline.bed.gz').touch()
            run(['bash','split_scripts/run-stage.sh','CaVEMan_flag','1'])
            self.assertTrue((work/'CaVEMan_flag.stage.tar.gz').is_file())
            tracker=json.loads((base/'tracker.json').read_text())
            self.assertEqual(tracker['active'],0)
            self.assertEqual(tracker['peak'],workers)
            # Context remained immutable throughout all successful and failed tasks.
            for name,digest in plan['seed_sha256'].items():
                self.assertEqual(scatter.sha(store/'context'/name),digest)


if __name__=='__main__': unittest.main(verbosity=2)
