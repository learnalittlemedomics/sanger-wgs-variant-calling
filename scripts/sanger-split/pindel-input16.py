#!/usr/bin/env python3
"""Official cgpPindel stages, isolated per chromosome; compatible with Python 3.5."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

PREFIX = 'Sanger_CGP_Pindel_Implement_'
OOM = re.compile(r'out.of.memory|oom[_ -]?kill|cannot allocate memory|std::bad_alloc|MemoryError|killed by signal\s*9', re.I)

def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, sort_keys=True, indent=2)+'\n')

def read(path):
    return json.loads(path.read_text())

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def cp(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        raise RuntimeError('Duplicate output: '+str(dst))
    shutil.copy2(str(src), str(dst))

def clean(path):
    if path.is_symlink():
        raise RuntimeError('Refusing to remove staged input: '+str(path))
    if path.exists():
        shutil.rmtree(str(path))

def need(path):
    if not path.is_file():
        raise RuntimeError('Missing output/marker: '+str(path))

def run(cmd, stage, index, threads, tmp):
    args = list(cmd)
    args[args.index('-c')+1] = str(threads)
    args += ['-p',stage,'-d']
    if index is not None:
        args += ['-i',str(index)]
    log = Path('timings')/('pindel.{}.{}.log'.format(stage,index or 0))
    log.parent.mkdir(exist_ok=True)
    dump(log.with_suffix('.command.json'), args)
    oom = False
    with log.open('w') as f:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in io.TextIOWrapper(proc.stdout, encoding='utf-8', errors='replace'):
            f.write(line); f.flush()
            print(line, end='', flush=True)
            oom = oom or bool(OOM.search(line))
        code = proc.wait()
    if code:
        for p in (tmp/'logs').glob('*'):
            if p.is_file():
                with p.open(errors='replace') as f:
                    oom = oom or any(OOM.search(line) for line in f)
        raise SystemExit(137 if oom or code in (-9,137,247) else (code if code > 0 else 1))

def options(cmd, tmp):
    return dict(tumour=cmd[cmd.index('-t')+1],normal=cmd[cmd.index('-n')+1],
                reference=cmd[cmd.index('-r')+1],exclude=cmd[cmd.index('-e')+1],tmp=str(tmp))

def layout(opts, full=False):
    opts = dict(opts, layout=full)
    output = subprocess.check_output(['perl','split_scripts/pindel-layout.pl'],
               input=json.dumps(opts).encode('utf-8'))
    return json.loads(output.decode('utf-8'))

def stage_inputs(paths, tmp):
    if len(paths) != 2:
        raise RuntimeError('Exactly two input parts required')
    reports = [read(p/'report.json') for p in paths]
    if sorted(r['index'] for r in reports) != [1,2]:
        raise RuntimeError('Duplicate/missing input role')
    versions = set(r['version'] for r in reports)
    if len(versions) != 1:
        raise RuntimeError('Mixed Pindel versions')
    for p,r in zip(paths,reports):
        if r.get('status') != 'PASS' or r.get('mode') != 'input':
            raise RuntimeError('Invalid input report')
        sample = r['sample_dir']
        if Path(sample).name != sample or sample in ('.','..'):
            raise RuntimeError('Unsafe sample directory')
        target = tmp/sample
        if target.exists() or target.is_symlink():
            raise RuntimeError('Tumour/normal sanitised sample names collide')
        target.symlink_to(str((p/'data').resolve()))
    return reports

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('mode',choices=['input','plan','chromosome','gather'])
    ap.add_argument('--command',required=True)
    ap.add_argument('--index',type=int,default=0)
    ap.add_argument('--inputs',type=Path,nargs='*',default=[])
    ap.add_argument('--plan',type=Path)
    ap.add_argument('--parts',type=Path,nargs='*',default=[])
    a=ap.parse_args()
    cmd=shlex.split(a.command)
    out=Path(cmd[cmd.index('-o')+1]);relative=out.relative_to(Path.cwd())
    if len(relative.parts)!=2 or relative.name!='pindel':
        raise RuntimeError('Unexpected Pindel output path')
    clean(out); out.mkdir(parents=True)
    tmp=out/'tmpPindel';tmp.mkdir()
    (tmp/'progress').mkdir(); (tmp/'logs').mkdir()
    opts=options(cmd,tmp)
    meta=layout(opts)
    if a.mode=='input':
        if a.index not in (1,2):raise RuntimeError('Input index must be 1 or 2')
        run(cmd,'input',a.index,16,tmp)
        need(tmp/'progress'/(PREFIX+'input.'+str(a.index)))
        sample=meta['tumour_dir' if a.index==1 else 'normal_dir']
        src=tmp/sample
        if not src.is_dir():raise RuntimeError('Input sample directory missing')
        result=Path('pindel_input_'+str(a.index));clean(result);result.mkdir()
        shutil.move(str(src),str(result/'data'))
        dump(result/'report.json',dict(mode='input',status='PASS',index=a.index,
             sample_dir=sample,version=meta['version'],tumour=meta['tumour_name'],normal=meta['normal_name']))
        return
    if a.mode in ('plan','chromosome'):
        reports=stage_inputs(a.inputs,tmp)
        for r in reports:
            if (r['version'],r['tumour'],r['normal']) != (meta['version'],meta['tumour_name'],meta['normal_name']):
                raise RuntimeError('Input sample/version mismatch')
        current=layout(opts,True)
        if current['jobs']<1 or current['jobs']!=len(current['chromosomes']):
            raise RuntimeError('No valid official chromosome jobs')
        if a.mode=='plan':
            current['inputs']=[dict(report=read(p/'report.json'),sha256=digest(p/'report.json')) for p in a.inputs]
            current['command']=cmd
            dump(Path('pindel-plan.json'),current)
            print('Pindel plan: {} chromosome jobs'.format(current['jobs']))
            return
        plan=read(a.plan)
        if current!=dict((k,v) for k,v in plan.items() if k not in ('inputs','command')):
            raise RuntimeError('Official chromosome mapping changed')
        if sorted(digest(p/'report.json') for p in a.inputs)!=sorted(x['sha256'] for x in plan['inputs']):
            raise RuntimeError('Input provenance mismatch')
        if not 1<=a.index<=plan['jobs']:raise RuntimeError('Chromosome index out of range')
        chrom=plan['chromosomes'][a.index-1]
        # Each task has all input filenames, so official sorted index mapping is unchanged.
        run(cmd,'pindel',a.index,1,tmp)
        need(tmp/'progress'/(PREFIX+'pindel.'+str(a.index)))
        for ext in ('D','SI'):need(tmp/'pout'/(chrom+'_'+chrom+'_'+ext))
        run(cmd,'pin2vcf',a.index,1,tmp)
        need(tmp/'progress'/(PREFIX+'pindel_to_vcf.'+str(a.index)))
        need(tmp/'vcf'/(chrom+'_pindel.vcf'))
        result=Path('pindel_part_'+str(a.index));clean(result);result.mkdir()
        shutil.move(str(tmp/'vcf'),str(result/'vcf'))
        dump(result/'report.json',dict(mode='chromosome',status='PASS',index=a.index,
             chromosome=chrom,plan_sha256=digest(a.plan)))
        return
    plan=read(a.plan);seen=set()
    if meta['version']!=plan['version']:raise RuntimeError('Pindel version changed')
    for part in a.parts:
        r=read(part/'report.json');i=r['index']
        if r.get('status')!='PASS' or r.get('mode')!='chromosome' or r.get('plan_sha256')!=digest(a.plan):
            raise RuntimeError('Failed/mixed chromosome result')
        if i in seen or not 1<=i<=plan['jobs'] or r['chromosome']!=plan['chromosomes'][i-1]:
            raise RuntimeError('Duplicate/invalid chromosome result')
        seen.add(i)
        need(part/'vcf'/(r['chromosome']+'_pindel.vcf'))
        for p in (part/'vcf').iterdir():
            if not p.is_file():raise RuntimeError('Unexpected chromosome output')
            cp(p,tmp/'vcf'/p.name)
    if seen!=set(range(1,plan['jobs']+1)):raise RuntimeError('Incomplete chromosome coverage')
    run(cmd,'merge',None,1,tmp)
    need(tmp/'progress'/(PREFIX+'merge_and_bam.0'))
    run(cmd,'flag',None,1,tmp)
    need(tmp/'progress'/(PREFIX+'flag.0'))
    stub=meta['tumour_name']+'_vs_'+meta['normal_name']
    for ext in ('.flagged.vcf.gz','.flagged.vcf.gz.tbi','.germline.bed'):need(out/(stub+ext))
    # Original downstream contract: retain flagged VCF, alignment outputs, gzip/index germline BED.
    bed=out/(stub+'.germline.bed');sorted_bed=out/(stub+'.germline.sorted.bed')
    with sorted_bed.open('wb') as f:
        subprocess.run(['sort','-k1,1','-k2,2n','-k3,3n',str(bed)],stdout=f,check=True)
    with (out/(stub+'.germline.bed.gz')).open('wb') as f:
        subprocess.run(['bgzip','-c',str(sorted_bed)],stdout=f,check=True)
    subprocess.run(['tabix','-p','bed',str(out/(stub+'.germline.bed.gz'))],check=True)
    bed.unlink();sorted_bed.unlink()
    # Same retention as official cleanup; don't call private Perl internals.
    shutil.move(str(tmp/'logs'),str(out/'logs'))
    clean(tmp)
    for p in out.glob('*.vcf.gz*'):
        if not p.name.endswith(('.flagged.vcf.gz','.flagged.vcf.gz.tbi')):p.unlink()
    dump(out/'scatter-report.json',dict(status='PASS',chromosomes=plan['chromosomes'],plan_sha256=digest(a.plan)))
    subprocess.run(['tar','-czf','cgpPindel.stage.tar.gz',str(relative),'timings'],check=True)

if __name__=='__main__':main()
