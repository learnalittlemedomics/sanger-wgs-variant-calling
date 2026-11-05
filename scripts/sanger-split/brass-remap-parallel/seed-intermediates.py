#!/usr/bin/env python3
"""Copy an explicitly verified, inactive BRASS directory into a new workdir."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil

MARKER = 'tmpBrass/progress/Sanger_CGP_Brass_Implement_filter.abs_bkp'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def seed(source, destination, r4_name, expected):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError('Source and destination must be separate directory trees')
    if Path(r4_name).name != r4_name or not r4_name.endswith('.r4'):
        raise ValueError('r4-name must be a plain .r4 filename')
    if not re.match(r'^[0-9a-f]{64}$', expected):
        raise ValueError('Expected SHA256 must be 64 lowercase hexadecimal characters')
    def verify(directory):
        if not (directory / MARKER).is_file() or digest(directory / r4_name) != expected:
            raise ValueError('Missing completion marker or r4 checksum mismatch')
    verify(source)
    if destination.exists():
        verify(destination)
        return {'status': 'REUSED', 'r4_sha256': expected}
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = destination.with_name(destination.name + '.seed-' + str(os.getpid()))
    if staging.exists():
        raise ValueError('Staging directory already exists; inspect it before retrying')
    shutil.copytree(str(source), str(staging), symlinks=True)
    verify(staging)
    os.rename(str(staging), str(destination))
    return {'status': 'COPIED', 'r4_sha256': expected}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--destination', required=True)
    parser.add_argument('--r4-name', required=True)
    parser.add_argument('--expected-r4-sha256', required=True)
    args = parser.parse_args()
    print(json.dumps(seed(args.source, args.destination, args.r4_name, args.expected_r4_sha256)))

if __name__ == '__main__':
    main()
