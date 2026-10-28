#!/usr/bin/env python3
"""Schedule the installed BRASS implementation without changing calling thresholds.
Compatible with the production container's Python 3.5.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tarfile


def digest_stream(stream):
    h = hashlib.sha256()
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        h.update(block)
    return h.digest()


def extract(archive):
    # Parallel branches have separate writable directories; join only identical
    # duplicates. Never silently overwrite a different branch's file/progress.
    with tarfile.open(str(archive), 'r:*') as tf:
        for m in tf.getmembers():
            p = Path(m.name)
            if p.is_absolute() or '..' in p.parts or m.islnk() or not (m.isdir() or m.isfile() or m.issym()):
                raise RuntimeError('Unsafe archive member: ' + m.name)
            if any(parent.is_symlink() for parent in p.parents):
                raise RuntimeError('Archive writes through symlink: ' + m.name)
            if p.exists() or p.is_symlink():
                same = m.isdir() and p.is_dir() and not p.is_symlink()
                if m.issym():
                    same = p.is_symlink() and os.readlink(str(p)) == m.linkname
                elif m.isfile() and p.is_file() and not p.is_symlink():
                    with p.open('rb') as existing, tf.extractfile(m) as incoming:
                        same = digest_stream(existing) == digest_stream(incoming)
                if not same:
                    raise RuntimeError('Conflicting branch output: ' + m.name)
            else:
                tf.extract(m, '.')


def pack(name, paths, compressed=False):
    with tarfile.open(name, 'w:gz' if compressed else 'w', dereference=False) as tf:
        for path in paths:
            path = Path(path)
            if path.exists() or path.is_symlink():
                tf.add(str(path), arcname=str(path))


def patched_input_env(cpus):
    if cpus < 30:
        raise RuntimeError('BRASS input10 requires 30 allocated CPUs')
    module = subprocess.check_output(['perl', '-MSanger::CGP::Brass::Implement', '-e',
                                      'print $INC{"Sanger/CGP/Brass/Implement.pm"}']).decode()
    source = Path(module).read_text()
    line = '$helper_threads = 6 if($helper_threads > 6);'
    if source.count(line) != 1:
        raise RuntimeError('Installed BRASS helper cap differs from audited source; refusing patch')
    dst = Path('brass_input_lib/Sanger/CGP/Brass/Implement.pm')
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(source.replace(line, '$helper_threads = 10 if($helper_threads > 10);'))
    env = os.environ.copy()
    env['PERL5LIB'] = str(Path('brass_input_lib').resolve()) + (':' + env['PERL5LIB'] if env.get('PERL5LIB') else '')
    # Preload before brass.pl's own 'use lib' can change module precedence.
    env['PERL5OPT'] = (env.get('PERL5OPT', '') + ' -MSanger::CGP::Brass::Implement').strip()
    loaded = subprocess.check_output(['perl', '-MSanger::CGP::Brass::Implement', '-e',
                                     'print $INC{"Sanger/CGP/Brass/Implement.pm"}'], env=env).decode()
    if Path(loaded).resolve() != dst.resolve():
        raise RuntimeError('Patched input module was not loaded')
    print('BRASS_INPUT_HELPER_THREADS=10 per input; allocation=30', flush=True)
    return env


def run(args):
    base = shlex.split(args.command)
    brass = Path(args.brass_dir)
    tmp = brass / 'tmpBrass'
    timings = Path('timings')
    timings.mkdir(exist_ok=True)
    env = patched_input_env(args.cpus) if args.stage == 'input' else None

    def invoke(stage, index=None):
        cmd = list(base)
        # Only these stages need ASCAT, as checked by brass.pl setup().
        if stage not in ('normcn', 'filter', 'grass', 'tabix') and '-ss' in cmd:
            n = cmd.index('-ss')
            del cmd[n:n+2]
        if index is not None and stage == 'assemble':
            cmd[cmd.index('-c') + 1] = '1'
        cmd += ['-p', stage]
        if index is not None:
            cmd += ['-i', str(index)]
        log = timings / ('BRASS.' + stage + ('.' + str(index) if index is not None else '') + '.log')
        with log.open('w') as out:
            out.write('COMMAND: ' + ' '.join(shlex.quote(x) for x in cmd) + '\n')
            out.flush()
            timed = ['/usr/bin/time', '-v'] + cmd if Path('/usr/bin/time').is_file() else cmd
            code = subprocess.call(timed, stdout=out, stderr=subprocess.STDOUT, env=env)
        if code:
            print('BRASS subcommand failed; log: ' + str(log), file=sys.stderr)
            print(log.read_text(errors='replace')[-12000:], file=sys.stderr)
            tail = log.read_text(errors='replace')[-65536:]
            if code in (-9, 137, 140, 143, 247) or re.search(r'out.of.memory|oom[_ -]?kill|cannot allocate memory|MemoryError', tail, re.I):
                print('SCATTER_RESOURCE_FAILURE: ' + str(log), file=sys.stderr)
                raise SystemExit(137)
            raise subprocess.CalledProcessError(code, cmd)

    if args.stage in ('input', 'cover', 'group', 'isize'):
        if args.stage == 'input':
            if args.group not in (1, 2):
                raise RuntimeError('Input role must be 1 (tumour) or 2 (normal)')
            invoke('input', args.group)
        else:
            invoke(args.stage)
    elif args.stage == 'normcn':
        invoke('merge')
        invoke('normcn')
    elif args.stage == 'plan':
        invoke('filter')
        invoke('split')
        files = sorted((tmp / 'split').glob('split.*'))
        expected = ['split.{:07d}'.format(i) for i in range(len(files))]
        if not files or [p.name for p in files] != expected:
            raise RuntimeError('Missing/non-contiguous official assembly splits')
        # Preserve the official 40-record split, including its empty-data split.
        plan = {'splits': len(files), 'groups': [list(range(i, min(i+6, len(files)+1)))
                                              for i in range(1, len(files)+1, 6)]}
        Path('brass-plan.json').write_text(json.dumps(plan, sort_keys=True))
        pack('brass-seed.tar', [tmp / 'split', 'brass-plan.json'])
    elif args.stage == 'assemble':
        plan = json.loads(Path('brass-plan.json').read_text())
        if not 1 <= args.group <= len(plan['groups']) or not 1 <= args.workers <= 6:
            raise RuntimeError('Invalid assembly group/workers')
        indices = plan['groups'][args.group-1]
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(invoke, 'assemble', i) for i in indices]
            for future in concurrent.futures.as_completed(futures):
                future.result()
        for i in indices:
            if not (tmp / 'assemble' / ('bedpe.{:07d}'.format(i-1))).is_file():
                raise RuntimeError('Assembly output missing for index ' + str(i))
        receipt = Path('brass-group-{}.json'.format(args.group))
        receipt.write_text(json.dumps(indices))
        pack('brass-part-{}.tar'.format(args.group),
             [tmp / 'assemble', tmp / 'progress', tmp / 'logs', timings, receipt])
        return
    elif args.stage == 'gather':
        plan = json.loads(Path('brass-plan.json').read_text())
        received = []
        for group, indices in enumerate(plan['groups'], 1):
            path = Path('brass-group-{}.json'.format(group))
            if not path.is_file() or json.loads(path.read_text()) != indices:
                raise RuntimeError('Missing/wrong assembly group ' + str(group))
            received.extend(indices)
        if sorted(received) != list(range(1, plan['splits']+1)):
            raise RuntimeError('Duplicate/missing assembly indices')
        for i in received:
            if not (tmp / 'assemble' / ('bedpe.{:07d}'.format(i-1))).is_file():
                raise RuntimeError('Missing assembly output ' + str(i))
        actual = sorted(p.name for p in (tmp / 'assemble').glob('bedpe.*'))
        if actual != ['bedpe.{:07d}'.format(i-1) for i in received]:
            raise RuntimeError('Unexpected assembly outputs')
        invoke('grass')
        invoke('tabix')
        pack('BRASS.stage.tar.gz', [brass, timings], compressed=True)
        return
    else:
        raise RuntimeError('Unknown BRASS stage: ' + args.stage)
    paths = [brass, timings]
    if args.stage == 'plan':
        paths += ['brass-plan.json']
    suffix = args.stage + ('_' + str(args.group) if args.stage == 'input' else '')
    pack('BRASS_{}.stage.tar'.format(suffix), paths)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='action')
    ex = sub.add_parser('extract')
    ex.add_argument('archive')
    r = sub.add_parser('run')
    r.add_argument('--stage', required=True)
    r.add_argument('--command', required=True)
    r.add_argument('--brass-dir', required=True)
    r.add_argument('--cpus', type=int, required=True)
    r.add_argument('--group', type=int, default=0)
    r.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    if args.action == 'extract':
        extract(args.archive)
    elif args.action == 'run':
        run(args)
    else:
        parser.error('subcommand required')


if __name__ == '__main__':
    main()
