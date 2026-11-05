#!/usr/bin/env python3
"""Public synthetic orchestration checks; these do not validate variant calls."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

REPO = Path(__file__).resolve().parents[2]
RELEASE = REPO / 'scripts/sanger-split/brass-remap-parallel'
ARRAYS = ['low_end_remap_score', 'high_end_remap_score', 'low_end_footprint', 'high_end_footprint']

# A fake interpreter implements only the worker's serialization protocol.
# Synthetic strings and integers deliberately replace biological scoring.
FAKE_PERL = r'''#!/usr/bin/env python3
import json,os,pathlib,sys,time
args=sys.argv[2:]; opts={}; pos=[]; i=0
while i<len(args):
    if args[i].startswith('-'):
        key=args[i][1:]
        if key=='noclean':opts[key]=True;i+=1
        else:opts[key]=args[i+1];i+=2
    else:pos.append(args[i]);i+=1
arrays=['low_end_remap_score','high_end_remap_score','low_end_footprint','high_end_footprint']
if 'load_remap_state' in opts:
    state=json.loads(pathlib.Path(opts['load_remap_state']).read_text())
    result={k:state[k] for k in arrays+['rg_id_of_read']}
    result['reads']=''.join(pathlib.Path(p).read_text() for p in state['all_reads_files'])
    pathlib.Path(pos[3]).write_text(json.dumps(result,sort_keys=True))
    pathlib.Path(opts['scores_output_file']).write_text(','.join(map(str,state[arrays[0]])))
else:
    rows=[line.split('\t') for line in pathlib.Path(pos[0]).read_text().splitlines() if not line.startswith('#')]
    if any(row[6]==os.environ.get('STUB_FAIL_GROUP') for row in rows):sys.exit(17)
    time.sleep(float(os.environ.get('STUB_SLEEP','0')))
    state={k:[] for k in arrays};state['rg_id_of_read']={};reads=[]
    for row in rows:
        value=int(row[1]);group=row[6]
        for k in arrays:state[k].append(value)
        state['rg_id_of_read']['shared-read']=group
        state['rg_id_of_read']['read-'+str(value)]=group
        reads.append('>'+group+'\nACGT\n')
    fasta=pathlib.Path(opts['tmpdir'])/'synthetic.fa';fasta.write_text(''.join(reads))
    state['all_reads_files']=[str(fasta)]
    pathlib.Path(opts['remap_state']).write_text(json.dumps(state,sort_keys=True))
'''

class ReleaseScope(unittest.TestCase):
    def test_exact_deployed_release(self):
        record=json.loads((REPO/'docs/benchmarks/brass-remap-parallel-2026-11-06.json').read_text())
        for name, expected in record['published_release_sha256'].items():
            # Git may convert line endings on Windows; the deployed files use LF.
            actual=(RELEASE/name).read_bytes().replace(b'\r\n',b'\n')
            self.assertEqual(hashlib.sha256(actual).hexdigest(),expected,name)

    def test_worker_only_adds_state_boundary(self):
        text=(RELEASE/'worker.pl').read_text(encoding='utf-8')
        text=text.replace('use strict;\nuse JSON::PP;', 'use strict;', 1)
        text=text.replace('my $REMAP_STATE = q{};\nmy $LOAD_REMAP_STATE = q{};\n','',1)
        text=text.replace("  'remap_state=s' => \\$REMAP_STATE,\n  'load_remap_state=s' => \\$LOAD_REMAP_STATE,\n",'',1)
        start=text.index('if ($LOAD_REMAP_STATE ne q{}) {')
        loop=text.index('for my $i(0..$#regions) {',start)
        text=text[:start]+text[loop:]
        start=text.index('}\nif ($REMAP_STATE ne q{}) {')
        end=text.index('#die "REQUSESTED EXIT POINT";',start)
        text=text[:start]+text[end:]
        self.assertEqual(hashlib.sha256(text.encode()).hexdigest(),
                         'ff7e5d886547e2c52aded3c3f1d694139129f615a42ecdfbd9d2502142b99503')

@unittest.skipUnless(os.name=='posix','Controller uses Linux process groups and flock')
class ControllerBehavior(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name);self.bin=self.base/'bin';self.bin.mkdir()
        (self.bin/'perl').write_text(FAKE_PERL);(self.bin/'perl').chmod(0o755)
        self.worker=self.base/'worker.pl';self.worker.write_text('synthetic protocol stub\n')
        self.input=self.base/'input.r4';self.groups=self.base/'groups.tsv'
        self.make_rows(['g0','g1','g2','g0','g4','g5','g6','g7','g8'])
        for name in ['input.bam','ref.fa','viral.2bit']+['bacteria.'+str(i)+'.fa.2bit' for i in range(1,5)]:
            (self.base/name).write_text('synthetic\n')
        runner=self.base/'run.py'
        runner.write_text("import importlib.util,signal,sys\n"
            "spec=importlib.util.spec_from_file_location('controller',sys.argv[1]);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)\n"
            "m.WORKER=m.pathlib.Path(sys.argv[2])\n"
            "signal.signal(signal.SIGTERM,m.interrupted);signal.signal(signal.SIGINT,m.interrupted)\n"
            "try:m.main(sys.argv[3:])\n"
            "except BaseException:m.terminate_children();raise\n")
        self.cmd=[sys.executable,str(runner),str(RELEASE/'controller.py'),str(self.worker),
            '-groups_file',str(self.groups),'-scores_output_file',str(self.base/'result.scores'),
            '-tmpdir',str(self.base/'tmp'),'-virus_db',str(self.base/'viral.2bit'),
            '-bacterial_db_stub',str(self.base/'bacteria'),'-search_cores','2',
            str(self.input),str(self.base/'input.bam'),str(self.base/'ref.fa'),str(self.base/'result.r5')]
        self.env=dict(os.environ,PATH=str(self.bin)+os.pathsep+os.environ['PATH'],BRASS_REMAP_CHUNK_SIZE='2')
        self.root=self.base/'tmp/parallel-v1'

    def make_rows(self, groups):
        self.ids=groups
        rows=['chrSynthetic\t%d\t%d\tchrSynthetic\t1\t2\t%s\t3\t+\t-\n'%(i,i+1,g) for i,g in enumerate(groups)]
        self.input.write_text('# synthetic test\n'+''.join(rows))
        self.groups.write_text(''.join(reversed(rows)))

    def run_it(self, **env):
        return subprocess.run(self.cmd,env=dict(self.env,**env),stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30)

    def assert_success(self,p):
        self.assertEqual(p.returncode,0,p.stderr.decode(errors='replace'))

    def test_order_shared_reads_group_routing_and_cached_result(self):
        self.cmd[self.cmd.index('-search_cores')+1]='30'
        first=self.run_it(STUB_SLEEP='0.1');self.assert_success(first)
        self.assertIn(b'workers=15 threads=2',first.stdout)
        result=json.loads((self.base/'result.r5').read_text())
        for key in ARRAYS:self.assertEqual(result[key],list(range(len(self.ids))))
        self.assertEqual(result['rg_id_of_read']['shared-read'],self.ids[-1])
        self.assertEqual(result['reads'],''.join('>'+g+'\nACGT\n' for g in self.ids))
        all_groups=self.groups.read_text().splitlines(True)
        for chunk in self.root.glob('chunk-*'):
            needed={line.split('\t')[6] for line in (chunk/'input.r4').read_text().splitlines()}
            self.assertEqual((chunk/'groups.tsv').read_text(),''.join(line for line in all_groups if line.split('\t')[6] in needed))
        receipts={p:p.stat().st_mtime_ns for p in self.root.glob('chunk-*/complete.json')}
        p=self.run_it();self.assert_success(p);self.assertIn(b'COMPLETE_CACHED',p.stdout)
        self.assertTrue(all(p.stat().st_mtime_ns==mtime for p,mtime in receipts.items()))

    def test_actual_interrupt_then_resume_reuses_complete_chunks(self):
        p=subprocess.Popen(self.cmd,env=dict(self.env,STUB_SLEEP='0.25'),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            deadline=time.monotonic()+10
            while time.monotonic()<deadline:
                if list(self.root.glob('chunk-*/complete.json')):break
                if p.poll() is not None:self.fail('Exited before checkpoint')
                time.sleep(0.02)
            self.assertTrue(list(self.root.glob('chunk-*/complete.json')))
            p.terminate();p.communicate(timeout=10);self.assertNotEqual(p.returncode,0)
        finally:
            if p.poll() is None:p.kill();p.communicate()
        receipts={f:f.stat().st_mtime_ns for f in self.root.glob('chunk-*/complete.json')}
        self.assert_success(self.run_it())
        self.assertTrue(all(f.stat().st_mtime_ns==mtime for f,mtime in receipts.items()))
        self.assertEqual(json.loads((self.base/'result.r5').read_text())[ARRAYS[0]],list(range(len(self.ids))))

    def test_corrupted_checkpoint_is_rejected(self):
        self.assert_success(self.run_it());(self.root/'complete.json').unlink()
        state=next(self.root.glob('chunk-*/state.json'));state.write_bytes(state.read_bytes()+b' ')
        p=self.run_it();self.assertNotEqual(p.returncode,0);self.assertIn(b'Corrupt checkpoint',p.stderr)
        self.assertFalse((self.root/'complete.json').exists())

    def test_changed_input_signature_is_rejected(self):
        self.assert_success(self.run_it());self.input.write_text(self.input.read_text()+'# changed\n')
        p=self.run_it();self.assertNotEqual(p.returncode,0);self.assertIn(b'different inputs/code/options',p.stderr)

    def test_child_failure_never_marks_success(self):
        p=self.run_it(STUB_FAIL_GROUP='g2');self.assertNotEqual(p.returncode,0)
        self.assertFalse((self.root/'complete.json').exists());self.assertFalse((self.base/'result.r5').exists())
        self.assertEqual(json.loads((self.root/'progress.json').read_text())['stage'],'failed')

    def test_empty_input(self):
        self.make_rows([]);self.assert_success(self.run_it())
        result=json.loads((self.base/'result.r5').read_text())
        self.assertEqual(result[ARRAYS[0]],[]);self.assertEqual(result['reads'],'')

    def test_damaged_final_output_is_rejected(self):
        self.assert_success(self.run_it());(self.base/'result.scores').write_text('damaged')
        p=self.run_it();self.assertNotEqual(p.returncode,0);self.assertIn(b'Final checkpoint/output mismatch',p.stderr)

class SeedBehavior(unittest.TestCase):
    def setUp(self):
        spec=importlib.util.spec_from_file_location('seed',RELEASE/'seed-intermediates.py')
        self.module=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.module)
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name);self.source=self.base/'source';self.source.mkdir()
        marker=self.source/self.module.MARKER;marker.parent.mkdir(parents=True);marker.touch()
        (self.source/'synthetic.r4').write_bytes(b'synthetic\n');self.expected=hashlib.sha256(b'synthetic\n').hexdigest()
        self.dest=self.base/'destination'

    def test_copy_and_verify_existing_destination(self):
        self.assertEqual(self.module.seed(self.source,self.dest,'synthetic.r4',self.expected)['status'],'COPIED')
        self.assertEqual(self.module.seed(self.source,self.dest,'synthetic.r4',self.expected)['status'],'REUSED')
        (self.dest/'synthetic.r4').write_bytes(b'damaged')
        with self.assertRaises(ValueError):self.module.seed(self.source,self.dest,'synthetic.r4',self.expected)

    def test_bad_hash_missing_marker_and_overlapping_paths_rejected(self):
        with self.assertRaises(ValueError):self.module.seed(self.source,self.dest,'synthetic.r4','0'*64)
        with self.assertRaises(ValueError):self.module.seed(self.source,self.source/'nested','synthetic.r4',self.expected)
        (self.source/self.module.MARKER).unlink()
        with self.assertRaises(ValueError):self.module.seed(self.source,self.dest,'synthetic.r4',self.expected)
        self.assertFalse(self.dest.exists())

if __name__=='__main__':unittest.main(verbosity=2)
