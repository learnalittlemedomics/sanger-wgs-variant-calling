#!/usr/bin/env python3
"""Flag-only scheduling adapter; scientific flags remain in the container tools."""
import sys, os, shlex, shutil, subprocess, tarfile, json
from pathlib import Path
mode, command, ident, input_file = sys.argv[1:]
args = shlex.split(command)
opts = dict(zip(args[1::2], args[2::2]))
cave = Path(opts['-o']); tmp = cave/'tmpCaveman'

def pack(name, paths):
    with tarfile.open(name, 'w:gz') as out:
        for path in paths:
            out.add(str(path), arcname=str(path.relative_to(Path.cwd())) if path.is_absolute() else str(path))

def run(cmd):
    print('COMMAND: '+ ' '.join(shlex.quote(x) for x in cmd), flush=True)
    subprocess.check_call(cmd)

def flag_one(input_file):
    ident = Path(input_file).name.rsplit('.',1)[1]
    # Derive sample prefix from the official splitter output filename.
    name = Path(input_file).name
    suffix = '.muts.ids.vcfsplit.'+ident
    if not name.endswith(suffix): raise ValueError('Unexpected split filename: '+name)
    output = tmp/(name[:-len(suffix)]+'.flagged.muts.vcf.'+ident)
    cmd = ['perl', shutil.which('cgpFlagCaVEMan.pl'), '-i', input_file, '-o', str(output)]
    resolved = json.loads(Path('flag-options.json').read_text())
    for src,dst in [('species','-s'),('tumbam','-m'),('normbam','-n'),('flag-bed','-b'),('germindel','-g'),('unmatchedvcf','-umv'),('reference','-ref'),('seqType','-t'),('species-assembly','-sa'),('flagConfig','-c'),('flagToVcfConfig','-v'),('annot-bed','-ab'),('apid','-p')]:
        if resolved.get(src) is not None: cmd += [dst,str(resolved[src])]
    logs = Path('timings')/('flag.'+ident+'.log')
    with logs.open('w') as handle:
        handle.write('COMMAND: '+ ' '.join(shlex.quote(x) for x in cmd)+'\n'); handle.flush()
        subprocess.check_call(cmd, stdout=handle, stderr=subprocess.STDOUT)
    if not output.is_file(): raise RuntimeError('Flag output missing')
    def records(path):
        with path.open() as handle:
            return sum(1 for line in handle if line.strip() and not line.startswith('#'))
    if records(Path(input_file)) != records(output):
        raise RuntimeError('Flag output record count differs from input')
    marker = tmp/'progress'/('Sanger_CGP_Caveman_Implement_caveman_flag.'+ident)
    marker.parent.mkdir(parents=True, exist_ok=True); marker.touch()
    return [output, marker, logs]

if mode == 'worker':
    from concurrent.futures import ThreadPoolExecutor, as_completed
    ident = str(int(ident))
    if int(ident) < 1: raise ValueError('Invalid group ID')
    tmp.mkdir(parents=True, exist_ok=True)
    Path('timings').mkdir(exist_ok=True)
    with tarfile.open(input_file, 'r:gz') as archive:
        for m in archive.getmembers():
            if not m.isfile() or '/' in m.name or not m.name.rsplit('.',1)[-1].isdigit():
                raise RuntimeError('Invalid flag group member')
            target = Path('flag-input')/m.name
            target.parent.mkdir(exist_ok=True)
            with target.open('wb') as handle: shutil.copyfileobj(archive.extractfile(m),handle)
    inputs = sorted(Path('flag-input').glob('*.muts.ids.vcfsplit.*'), key=lambda f:int(f.name.rsplit('.',1)[1]))
    if not inputs: raise RuntimeError('Empty flag group')
    workers = int(os.environ.get('FLAG_WORKERS','30'))
    if workers < 1 or workers > 30: raise ValueError('Invalid flag worker count')
    outputs = []; failures = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(flag_one,str(f)) for f in inputs]
        for future in as_completed(futures):
            try: outputs.extend(future.result())
            except Exception as error:
                failures.append(error)
    if failures:
        for error in failures: print(str(error), file=sys.stderr)
        sys.exit(137 if any(isinstance(e,subprocess.CalledProcessError) and e.returncode in [-9,137] for e in failures) else 1)
    pack('flag-part-'+ident+'.tar.gz', outputs)

else:
    source = Path(shutil.which('caveman.pl')).read_text()
    old = "$threads->run($options->{'vcf_split_count'}, 'caveman_flag', $options);"
    if source.count(old) != 1: raise RuntimeError('Unsupported caveman.pl: flag dispatch differs')
    if mode == 'plan':
        import re
        source, n = re.subn(r'(const my \$SPLIT_LINE_COUNT\s*=>\s*)1000(\s*;)', r'\g<1>200\2', source)
        if n != 1: raise RuntimeError('Unsupported flag split constant')
        source = source.replace(old, "require JSON::PP; open(my $fh, '>', 'flag-options.json') or die $!; print $fh JSON::PP::encode_json({ map { $_ => $options->{$_} } qw(species tumbam normbam flag-bed germindel unmatchedvcf reference seqType species-assembly flagConfig flagToVcfConfig annot-bed apid) }); close($fh); exit 0;")
    elif mode == 'gather':
        inputs = list(tmp.glob('*.muts.ids.vcfsplit.*'))
        if not inputs: raise RuntimeError('Missing flag plan')
        expected = set()
        for f in inputs:
            i = f.name.rsplit('.',1)[1]; expected.add(i)
            output = tmp/(f.name.split('.muts.ids.vcfsplit.')[0]+'.flagged.muts.vcf.'+i)
            marker = tmp/'progress'/('Sanger_CGP_Caveman_Implement_caveman_flag.'+i)
            if not output.is_file() or not marker.is_file(): raise RuntimeError('Missing shard '+i)
        actual = set(f.name.rsplit('.',1)[1] for f in tmp.glob('*.flagged.muts.vcf.*') if f.name.rsplit('.',1)[1].isdigit())
        if actual != expected: raise RuntimeError('Unexpected flag shard set')
        source = source.replace(old, '# All shard receipts verified before official concat/zip/cleanup')
    else: raise ValueError(mode)
    local = Path('caveman-flag-entry.pl'); local.write_text(source)
    run(['perl',str(local)] + args[1:])
    if mode == 'plan':
        splits = sorted(tmp.glob('*.muts.ids.vcfsplit.*'))
        if not splits: raise RuntimeError('Official splitter produced no output')
        # Round-robin genomic shards across groups to spread costly regions.
        splits.sort(key=lambda f:int(f.name.rsplit('.',1)[1]))
        groups = (len(splits)+299)//300
        Path('flag-groups').mkdir()
        for g in range(groups):
            with tarfile.open('flag-groups/flag-group-'+str(g+1)+'.tar.gz','w:gz') as archive:
                for f in splits[g::groups]: archive.add(str(f),arcname=f.name)
        pack('flag-plan.tar.gz', [cave])
    else:
        pack('CaVEMan_flag.stage.tar.gz', [cave,Path('timings')])
