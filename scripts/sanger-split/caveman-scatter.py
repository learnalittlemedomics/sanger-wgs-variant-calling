#!/usr/bin/env python3
"""Production CaVEMan scatter/gather. Runs INSIDE the existing Sanger container.
The genomic splitList and scientific options are unchanged. Only scheduling changes.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

PREFIX = 'Sanger_CGP_Caveman_Implement_caveman_'
SUFFIXES = ('muts.vcf.gz', 'snps.vcf.gz', 'no_analysis.bed')
OOM = re.compile(r'out.of.memory|oom[_ -]?kill|cannot allocate memory|std::bad_alloc|'
                 r'MemoryError|killed by signal\s*9|signal\s*[=:]?\s*9\b|exit(?:ed)?(?:\s+with)?(?:\s+(?:code|status))?\s*[=:]?\s*137\b', re.I)


def read_rows(path):
    rows = []
    names = set()
    for n, line in enumerate(path.read_text().splitlines(), 1):
        chrom, start, end = line.split()
        if int(start) < 0 or int(end) <= int(start) or '/' in chrom or chrom in ('.', '..'):
            raise ValueError('Invalid splitList row: ' + line)
        stem = '{}_{}'.format(int(start)+1, int(end))
        if (chrom, stem) in names:
            raise ValueError('Duplicate genomic shard: ' + line)
        names.add((chrom, stem))
        rows.append((n, chrom, stem))
    if not rows:
        raise ValueError('Empty splitList')
    return rows


def group_count(total, minutes, target_hours):
    if not math.isfinite(minutes) or not math.isfinite(target_hours) or minutes <= 0 or target_hours <= 0:
        raise ValueError('Timing parameters must be finite positive values')
    return max(1, min(total, math.ceil(total * minutes / (60 * target_hours))))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def copy_file(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(source), str(target))


def prepare(args):
    root = Path.cwd()
    cave = Path(args.cave_dir)
    relative = cave.relative_to(root)
    cmd = shlex.split(args.command)
    cmd[cmd.index('-t')+1] = '1'
    # split emits per-contig files; the official wrapper must assemble splitList.
    run_command(cmd, 'split_concat', ['-i', '1'], root / 'split-concat.log')
    rows = read_rows(cave / 'tmpCaveman/splitList')
    ctx = root / 'scatter_context'
    if ctx.is_symlink():
        raise RuntimeError('Refusing to overwrite staged context')
    if ctx.exists():
        shutil.rmtree(str(ctx))  # This is an output owned by the current plan task.
    ctx.mkdir()
    seed = ctx / 'seed'
    for name in ('alg_bean', 'splitList', 'caveman.cfg.ini'):
        copy_file(cave / 'tmpCaveman' / name, seed / 'tmpCaveman' / name)
    progress = cave / 'tmpCaveman/progress'
    for path in progress.glob('*'):
        if path.is_file() and any(path.name.startswith(PREFIX+s) for s in ('setup', 'split')):
            copy_file(path, seed / 'tmpCaveman/progress' / path.name)
    for name in ('norm.cn.bed', 'tum.cn.bed'):
        copy_file(root / 'tmp' / name, ctx / 'cn' / name)
    links = {}
    for key, bam, index in (('-tb', args.tumour, args.tidx), ('-nb', args.normal, args.nidx)):
        destination = Path(cmd[cmd.index(key)+1]).relative_to(root)
        # The context tuple stages these exact basenames in each task directory.
        links[str(destination)] = Path(bam).name
        links[str(destination)+'.bas'] = Path(bam).name+'.bas'
        index_source = Path(index).name
        index_destination = str(destination)+('.crai' if str(destination).endswith('.cram') else '.bai')
        if index_source.endswith('.csi'):
            index_destination = str(destination)+'.csi'
        links[index_destination] = index_source
    if not 1 <= args.workers_per_job <= 30:
        raise ValueError('workers-per-job must be in 1..30')
    plan = dict(format_version=2, workers_per_job=args.workers_per_job, total=len(rows), cave_relative=str(relative), command=cmd,
                links=links, mstep_minutes=args.mstep_minutes, estep_minutes=args.estep_minutes,
                target_hours=args.target_hours,
                mstep_groups=group_count(len(rows), args.mstep_minutes, args.target_hours),
                estep_groups=group_count(len(rows), args.estep_minutes, args.target_hours),
                seed_sha256={p.relative_to(ctx).as_posix(): sha(p) for p in ctx.rglob('*') if p.is_file()})
    for mode in ('mstep', 'estep'):
        plan[mode+'_worker_groups'] = plan[mode+'_groups']
        plan[mode+'_groups'] = math.ceil(plan[mode+'_worker_groups'] / args.workers_per_job)
    write_json(ctx / 'plan.json', plan)
    write_json(root / 'scatter-plan.json', plan)
    print('SCATTER PLAN: {} shards, {} Mstep groups, {} Estep groups'.format(
        plan['total'], plan['mstep_groups'], plan['estep_groups']))


def run_command(cmd, stage, extra, log):
    args = cmd + ['-p', stage] + extra
    write_json(log.with_suffix('.command.json'), args)
    env = dict(os.environ, PCAP_THREADED_REM_LOGS='0')
    # Keep logs of failed children: the Perl wrapper can mask a killed child as exit 1.
    start = time.monotonic()
    oom_seen = False
    with log.open('w') as out:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                env=env)
        for line in io.TextIOWrapper(proc.stdout, encoding='utf-8', errors='replace'):
            out.write(line)
            print(line, end='', flush=True)
            oom_seen = oom_seen or bool(OOM.search(line))
        code = proc.wait()
    if code:
        # Include child logs, while ignoring empty generic "Errors from command" headings.
        for directory in (Path.cwd()/'logs', Path(cmd[cmd.index('-o')+1])/'logs',
                          Path(cmd[cmd.index('-o')+1])/'tmpCaveman/logs'):
            if directory.exists():
                for p in directory.rglob('*'):
                    if p.is_file():
                        with p.open(errors='replace') as child:
                            oom_seen = oom_seen or any(OOM.search(line) for line in child)
        if oom_seen or code in (-9, 137, 247):
            print('SCATTER_RESOURCE_FAILURE: memory/SIGKILL; request resource retry', file=sys.stderr)
            raise SystemExit(137)
        raise SystemExit(code if 0 < code < 256 else 1)
    return time.monotonic()-start


def validate_output(path):
    if not path.is_file():
        raise RuntimeError('Missing output: ' + str(path))
    if path.name.endswith(('.gz', '.covs')):
        with gzip.open(str(path), 'rb') as handle:
            for block in iter(lambda: handle.read(1048576), b''):
                pass


def collect_parts(parts, expected_mode, rows, tmp, plan_hash):
    expected = {n for n, _, _ in rows}
    seen = set()
    for part in parts:
        report = json.loads((part/'report.json').read_text())
        if report.get('mode') != expected_mode or report.get('status') != 'PASS' or report.get('plan_sha256') != plan_hash:
            raise RuntimeError('Wrong/failed upstream part or mixed sample plan: ' + str(part))
        indices = report['indices']
        if len(indices) != len(set(indices)) or seen.intersection(indices) or not set(indices) <= expected:
            raise RuntimeError('Duplicate/out-of-range shard indices: ' + str(part))
        seen.update(indices)
        selected = [r for r in rows if r[0] in set(indices)]
        suffixes = ('covs',) if expected_mode == 'mstep' else SUFFIXES
        for n, chrom, stem in selected:
            for suffix in suffixes:
                source = part/'results'/chrom/(stem+'.'+suffix)
                destination = tmp/'results'/chrom/source.name
                if destination.exists():
                    raise RuntimeError('Duplicate output: ' + str(destination))
                copy_file(source, destination)
            marker = PREFIX+expected_mode+'.'+str(n)
            if not (part/'progress'/marker).is_file():
                raise RuntimeError('Missing upstream completion marker: ' + marker)
            copy_file(part/'progress'/marker, tmp/'progress'/marker)
    if seen != expected:
        raise RuntimeError('Incomplete shard coverage: {}/{}'.format(len(seen), len(expected)))


def execute(args):
    root = Path.cwd()
    ctx = args.context
    plan_file = ctx/'plan.json'
    plan = json.loads(plan_file.read_text())
    plan_hash = sha(plan_file)
    for name, digest in plan['seed_sha256'].items():
        if sha(ctx/name) != digest:
            raise RuntimeError('Scatter context changed: ' + name)
    relative = Path(plan['cave_relative'])
    if relative.is_absolute() or '..' in relative.parts or len(relative.parts) != 2 or relative.name != 'caveman':
        raise RuntimeError('Unsafe cave_relative in plan')
    cave = root / relative
    if cave.is_symlink() or (root/'result').is_symlink():
        raise RuntimeError('Refusing to replace a staged input')
    # Nextflow can retry in the same directory. Restart only this task's outputs;
    # never reuse partial child markers or write into staged upstream directories.
    for owned in (cave, root/'result'):
        if owned.exists():
            shutil.rmtree(str(owned))
    shutil.copytree(str(ctx/'seed'), str(cave))
    tmp = cave/'tmpCaveman'
    (tmp/'results').mkdir(exist_ok=True)
    (tmp/'progress').mkdir(exist_ok=True)
    (tmp/'logs').mkdir(exist_ok=True)
    for name in ('norm.cn.bed', 'tum.cn.bed'):
        copy_file(ctx/'cn'/name, root/'tmp'/name)
    for destination, source in plan['links'].items():
        p = root/destination
        if not (root/source).is_file():
            raise RuntimeError('Missing staged input: ' + source)
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.is_symlink():
            p.unlink()
        elif p.exists():
            raise RuntimeError('Unexpected existing input link destination: ' + str(p))
        p.symlink_to(str(root/source))
    os.environ['REF_PATH'] = str(root/'reference_files/genome.fa')
    os.environ['REF_CACHE'] = str(root/'ref_cache/%2s/%2s/%s')
    version = subprocess.check_output(['perl', '-MSanger::CGP::Caveman', '-e', 'print $Sanger::CGP::Caveman::VERSION'], universal_newlines=True).strip()
    if version != '1.13.14':
        raise RuntimeError('Validated wrapper version is 1.13.14, got ' + version)
    rows = read_rows(tmp/'splitList')
    if len(rows) != plan['total']:
        raise RuntimeError('splitList count differs from plan')
    cmd = plan['command']
    result = root/'result'
    result.mkdir()
    timings = root/'timings'
    timings.mkdir(exist_ok=True)
    log = timings/('scatter.{}.{}.log'.format(args.mode, args.index))
    elapsed = 0
    if args.mode in ('estep', 'gather'):
        model = json.loads((args.model/'report.json').read_text())
        if model.get('mode') != 'merge' or model.get('status') != 'PASS' or model.get('plan_sha256') != plan_hash:
            raise RuntimeError('Wrong model for this plan')
        for name in ('cov_arr', 'prob_arr'):
            copy_file(args.model/name, tmp/name)
        for p in (args.model/'progress').glob('*'):
            copy_file(p, tmp/'progress'/p.name)
    if args.mode in ('mstep', 'estep'):
        groups = plan[args.mode+'_groups']
        if not 1 <= args.index <= groups:
            raise ValueError('Group index out of range')
        worker_groups = plan.get(args.mode+'_worker_groups', groups)
        packing = plan.get('workers_per_job', 1)
        lane_ids = list(range((args.index-1)*packing+1, min(args.index*packing, worker_groups)+1))
        selected = [r for r in rows if (r[0]-1) % worker_groups + 1 in lane_ids]
        if not 1 <= args.workers <= 30:
            raise ValueError('workers must be in 1..30')
        # All subprocesses read the same immutable model/config; each lane owns
        # disjoint genomic output names, progress markers and child logs.
        if any(str(x).endswith('.gz.fai') for x in cmd):
            raise RuntimeError('Parallel worker requires an uncompressed reference .fai')
        for _, chrom, _ in selected:
            (tmp/'results'/chrom).mkdir(exist_ok=True)
        started = time.monotonic()
        failures = []
        lane_timings = {}
        with ThreadPoolExecutor(max_workers=min(args.workers, len(lane_ids))) as pool:
            futures = {pool.submit(run_command, cmd, args.mode,
                       ['-i', str(i), '-l', str(worker_groups)],
                       timings/('scatter.{}.lane_{}.log'.format(args.mode, i))): i for i in lane_ids}
            for future in as_completed(futures):
                i = futures[future]
                try:
                    lane_timings[str(i)] = future.result()
                except SystemExit as exc:
                    failures.append(int(exc.code))
                except Exception as exc:
                    print('Lane {} failed: {}'.format(i, exc), file=sys.stderr)
                    failures.append(1)
        elapsed = time.monotonic()-started
        write_json(timings/'pool.json', dict(workers=min(args.workers, len(lane_ids)),
                   lane_ids=lane_ids, wall_seconds=elapsed, lane_seconds=lane_timings,
                   failures=failures))
        if failures:
            # Do not publish a partially successful pool. All children have ended.
            raise SystemExit(137 if all(c in (137, 140, 143, 247) for c in failures) else 1)
        for n, chrom, stem in selected:
            marker = tmp/'progress'/(PREFIX+args.mode+'.'+str(n))
            if not marker.is_file():
                raise RuntimeError('Missing completion marker: ' + str(marker))
            for ext in (('covs',) if args.mode == 'mstep' else SUFFIXES):
                validate_output(tmp/'results'/chrom/(stem+'.'+ext))
            copy_file(marker, result/'progress'/marker.name)
        shutil.move(str(tmp/'results'), str(result/'results'))
        # Models/configs remain in the task tree, but are not emitted per group.
        indices = [n for n, _, _ in selected]
    else:
        upstream_mode = 'mstep' if args.mode == 'merge' else 'estep'
        collect_parts(args.parts, upstream_mode, rows, tmp, plan_hash)
        if args.mode == 'merge':
            elapsed = run_command(cmd, 'merge', [], log)
            for name in ('cov_arr', 'prob_arr'):
                if not (tmp/name).is_file() or (tmp/name).stat().st_size == 0:
                    raise RuntimeError('Missing/empty model: ' + name)
                copy_file(tmp/name, result/name)
            shutil.copytree(str(tmp/'progress'), str(result/'progress'))
        else:
            elapsed = run_command(cmd, 'merge_results', [], log)
            elapsed += run_command(cmd, 'add_ids', [], timings/'scatter.add_ids.log')
            # Preserve the original CaVEMan archive contract for flag/annotation/repack.
            write_json(cave/'scatter-plan.json', plan)
            for part in args.parts:
                report = json.loads((part/'report.json').read_text())
                copy_file(part/'report.json', cave/'scatter-reports'/('estep.{}.json'.format(report['index'])))
                if (part/'pool.json').is_file():
                    copy_file(part/'pool.json', cave/'scatter-reports'/('estep.{}.pool.json'.format(report['index'])))
            subprocess.run(['tar', '-czf', 'CaVEMan.stage.tar.gz', plan['cave_relative'], 'timings'], check=True)
        indices = [n for n, _, _ in rows]
    if args.mode in ('mstep', 'estep'):
        copy_file(timings/'pool.json', result/'pool.json')
    write_json(result/'report.json', dict(mode=args.mode, index=args.index, indices=indices,
               total=len(rows), status='PASS', plan_sha256=plan_hash, wall_seconds=elapsed))
    if args.mode in ('mstep', 'estep', 'merge'):
        shutil.rmtree(str(cave))  # Emitted result owns everything needed downstream.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    p = sub.add_parser('prepare')
    p.add_argument('--command', required=True)
    p.add_argument('--cave-dir', required=True)
    for name in ('tumour', 'normal', 'tidx', 'nidx'):
        p.add_argument('--'+name, required=True)
    p.add_argument('--workers-per-job', type=int, default=6)
    p.add_argument('--target-hours', type=float, default=1.5)
    p.add_argument('--mstep-minutes', type=float, default=3.328)
    p.add_argument('--estep-minutes', type=float, default=15.83)
    p = sub.add_parser('run')
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--mode', required=True, choices=['mstep', 'merge', 'estep', 'gather'])
    p.add_argument('--context', type=Path, default=Path('scatter_context'))
    p.add_argument('--index', type=int, default=1)
    p.add_argument('--model', type=Path, default=Path('model'))
    p.add_argument('--parts', type=Path, nargs='*', default=[])
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare(args)
    elif args.action == 'run':
        execute(args)
    else:
        parser.error('choose prepare or run')


if __name__ == '__main__':
    main()
