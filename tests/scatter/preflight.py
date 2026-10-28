#!/usr/bin/env python3
"""Validate production Nextflow processes locally with stubs, including an OOM retry."""
import csv
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parent
if not shutil.which('nextflow'):
    sys.exit('nextflow 未找到：请先激活运行 Sanger 的 Nextflow 环境。')
with tempfile.TemporaryDirectory(prefix='sanger-scatter-preflight-') as tmp:
    cmd = ['nextflow', '-C', str(root/'smoke.config'), 'run', str(root/'smoke.nf'),
           '-stub-run', '-work-dir', str(Path(tmp)/'work')]
    result = subprocess.run(cmd, cwd=tmp)
    if result.returncode:
        print('Nextflow 本地预检未通过；生产代码尚未安装。', file=sys.stderr)
        log = Path(tmp)/'.nextflow.log'
        if log.exists():
            print('\n'.join(log.read_text(errors='replace').splitlines()[-80:]), file=sys.stderr)
        sys.exit(result.returncode)
    rows = list(csv.DictReader((Path(tmp)/'smoke-trace.tsv').open(), delimiter='\t'))
    good = [r for r in rows if r['status']=='COMPLETED']
    assert len(good)==8, ('Expected 1 plan + 2 Mstep + 1 merge + 3 Estep + 1 gather', good)
    retried = [r for r in good if 'CaVEMan_Mstep' in r['name'] and 'group_1' in r['name']]
    assert len(retried)==1 and int(retried[0]['attempt'])==4 and int(retried[0]['cpus'])==30, retried
    other = [r for r in good if 'CaVEMan_Estep' in r['name']]
    assert len(other)==3 and all(int(r['cpus'])==24 for r in other), other
    resource_files = list((Path(tmp)/'work').glob('*/*/result/resources.txt'))
    resources = [p.read_text() for p in resource_files]
    assert any('attempt=4 cpus=30' in t and 'workers=2' in t for t in resources), resources
    print('PASS: Nextflow graph, grouped staging, model/gather dependencies, 137 retries through 24/30/30/30 CPUs; final workers=2.')

# No real Slurm calls: all scheduler commands are shadowed by isolated local stubs.
import os
with tempfile.TemporaryDirectory(prefix='sanger-submit-preflight-') as tmp:
    tmp=Path(tmp)
    fake=tmp/'fake_slurm'
    shutil.copytree(root/'fake_slurm',fake)
    for p in fake.iterdir():p.chmod(0o755)
    state=tmp/'state'
    env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH'],SCATTER_FAKE_SLURM_STATE=str(state))
    cmd=['nextflow','-C',str(root/'submit-smoke.config'),'run',str(root/'submit-smoke.nf'),
         '-work-dir',str(tmp/'work')]
    result=subprocess.run(cmd,cwd=str(tmp),env=env)
    if result.returncode:
        log=tmp/'.nextflow.log'
        if log.exists():print('\n'.join(log.read_text(errors='replace').splitlines()[-80:]),file=sys.stderr)
        sys.exit('提交重试预检失败；生产代码尚未安装。请保留上面的错误信息。')
    assert int((state/'submissions').read_text())==5,'Expected 4 rejected submissions then one acceptance'
    rows=list(csv.DictReader((tmp/'submit-trace.tsv').open(),delimiter='\t'))
    good=[r for r in rows if r['status']=='COMPLETED']
    assert len(good)==1 and int(good[0]['attempt'])==1,good
    print('PASS: 4 quota rejections retried; fifth submission accepted; task attempt remains 1 (no OOM budget/resource escalation).')

with tempfile.TemporaryDirectory(prefix='pindel-preflight-') as tmp:
    result=subprocess.run(['nextflow','-C',str(root/'smoke.config'),'run',str(root/'pindel-smoke.nf'),'-stub-run','-work-dir',str(Path(tmp)/'work')],cwd=tmp)
    if result.returncode: sys.exit('Pindel Nextflow preflight failed; not installed.')
    rows=list(csv.DictReader((Path(tmp)/'smoke-trace.tsv').open(),delimiter='\t'))
    good=[r for r in rows if r['status']=='COMPLETED']
    assert len(good)==6,good
    retried=[r for r in good if 'PindelChromosome (chr1)' in r['name']]
    assert len(retried)==1 and int(retried[0]['attempt'])==2 and int(retried[0]['cpus'])==12,retried
    print('PASS: Pindel input/plan/chromosome/gather graph, staged directories, chromosome retry 8 -> 12 CPUs.')

with tempfile.TemporaryDirectory(prefix='dependency-preflight-') as tmp:
    work=Path(tmp)/'work'
    result=subprocess.run(['nextflow','-C',str(root/'smoke.config'),'run',str(root/'dependency-smoke.nf'),'-stub-run','-work-dir',str(work)],cwd=tmp)
    if result.returncode: sys.exit('Dependency graph preflight failed; not installed.')
    rows=list(csv.DictReader((Path(tmp)/'smoke-trace.tsv').open(),delimiter='\t'))
    good=[r for r in rows if r['status']=='COMPLETED']
    assert len(good)==11,good
    starts={}; ends={}
    for script in work.glob('*/*/.command.sh'):
        text=script.read_text(); d=script.parent
        key=('bas2' if 'sleep 8' in text else 'bas1' if 'sleep 1' in text else
             'input1' if 'mkdir -p pindel_input_1' in text else
             'input2' if 'mkdir -p pindel_input_2' in text else
             'prepare' if 'touch reference_files/genome.fa' in text else
             'cache' if 'touch ref_cache/mock' in text else None)
        if key:
            starts[key]=(d/'.command.begin').stat().st_mtime
            ends[key]=(d/'.exitcode').stat().st_mtime
    assert ends['bas1'] <= starts['input1'] < ends['bas2'], (starts,ends)
    assert ends['bas2'] <= starts['input2'], (starts,ends)
    assert starts['prepare'] < ends['bas2'] and starts['cache'] < ends['bas2'], (starts,ends)
    print('PASS: preparation/cache overlap BAS; Input1 starts before BAS2 finishes; each input waits for its own BAS; downstream context reused for all chromosomes.')

with tempfile.TemporaryDirectory(prefix='brass-preflight-') as tmp:
    result=subprocess.run(['nextflow','-C',str(root/'smoke.config'),'run',str(root/'brass-smoke.nf'),'-stub-run','-work-dir',str(Path(tmp)/'work')],cwd=tmp)
    if result.returncode: sys.exit('BRASS Nextflow preflight failed; not installed.')
    rows=list(csv.DictReader((Path(tmp)/'smoke-trace.tsv').open(),delimiter='\t'))
    good=[r for r in rows if r['status']=='COMPLETED']
    assert len(good)==10, good
    inputs=[r for r in good if 'BRASSInput30' in r['name']]
    assert len(inputs)==2 and all(int(r['cpus'])==30 for r in inputs), inputs
    assert len([r for r in good if 'BRASSAssemble' in r['name']])==2, good
    print('PASS: BRASS two input roles, independent coverage/isize, plan, two assembly groups and gather.')
