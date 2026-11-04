#!/usr/bin/env python3
"""Run real-data differential tests inside a compatible BRASS environment.

Supply a bounded .r3 fixture, not the unbounded whole-data quadratic oracle.
This tool never overwrites the source scripts or input alignments.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('original','optimized','alignment','reference','candidates','output'):
        ap.add_argument('--'+name, required=True, type=Path)
    ap.add_argument('--repeats', type=int, default=2)
    args=ap.parse_args()
    if args.repeats < 1: ap.error('--repeats must be positive')
    args.output.mkdir(parents=True, exist_ok=False)
    args.output = args.output.resolve()
    runs=[]
    for repeat in range(args.repeats):
        order=('original','optimized') if repeat%2==0 else ('optimized','original')
        for variant in order:
            script=getattr(args,variant).resolve()
            stem=args.output/('%s-%d'%(variant,repeat))
            result=stem.with_suffix('.r4')
            started=time.monotonic(); markers={}
            with stem.with_suffix('.stderr').open('w') as log:
                proc=subprocess.Popen(['perl',str(script),'-fasta',str(args.reference.resolve()),
                    '-out',str(result),str(args.alignment.resolve()),str(args.candidates.resolve())],
                    stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,universal_newlines=True)
                for line in proc.stderr:
                    log.write(line)
                    for marker in ('Collecting supporting','Resolving overlapping','Outputting read position'):
                        if marker in line: markers[marker]=time.monotonic()-started
                rc=proc.wait()
            stderr=stem.with_suffix('.stderr').read_text()
            # Strip only source-location suffixes; keep all warning text,
            # occurrence counts, read-fetch events and their order.
            stderr=re.sub(r' at '+re.escape(str(script))+r' line \d+\.', ' at <script> line <n>.', stderr)
            record={'variant':variant,'repeat':repeat,'exit_code':rc,'elapsed_seconds':time.monotonic()-started,
                    'sha256':hashlib.sha256(result.read_bytes()).hexdigest() if result.exists() else None,
                    'stderr_sha256':hashlib.sha256(stderr.encode()).hexdigest(),'markers':markers}
            if 'Outputting read position' in markers:
                record['overlap_seconds']=markers['Outputting read position']-markers['Resolving overlapping']
            runs.append(record)
            (args.output/'runs.json').write_text(json.dumps(runs,indent=2))
    passed=all(r['exit_code']==0 and r['sha256'] is not None for r in runs) and len({r['sha256'] for r in runs})==1 and len({r['stderr_sha256'] for r in runs})==1
    (args.output/'summary.json').write_text(json.dumps({'passed':passed,'runs':runs},indent=2))
    print('PASS' if passed else 'FAIL')
    return 0 if passed else 1

if __name__=='__main__': raise SystemExit(main())
