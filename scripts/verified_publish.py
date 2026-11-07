#!/usr/bin/env python3
"""Synchronous publication for local ICGC workflows (Python 3, POSIX)."""
import argparse
import base64
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import uuid


def stamp(path):
    s = path.stat()
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)


def digest(path):
    h = hashlib.md5()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def sync_dir(path):
    fd = os.open(str(path), os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def save(path, data):
    temp = path.with_name('.' + path.name + '.' + uuid.uuid4().hex)
    try:
        with temp.open('x') as stream:
            json.dump(data, stream, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(str(temp), str(path))
        sync_dir(path.parent)
    finally:
        if temp.exists():
            temp.unlink()


def safe_name(value):
    if not isinstance(value, str) or not value or value in ('.', '..') or '/' in value or '\\' in value:
        raise ValueError('Unsafe output name: ' + repr(value))
    return value


def plan(manifest):
    root = Path(manifest['work_root']).resolve(strict=True)
    output = Path(manifest['output']).resolve()
    if output == root or root in output.parents:
        raise ValueError('Published results must be outside the work directory')
    expected = manifest['expected']
    groups = manifest['groups']
    if not groups or len(groups) != len(expected) or sorted(g['role'] for g in groups) != sorted(expected):
        raise ValueError('Missing or duplicate publication groups')
    rows = []
    seen = set()
    for group in groups:
        role = safe_name(group['role'])
        if 'raw_files' in group:
            if not group['raw_files']:
                raise ValueError('Empty supplementary output group')
            for entry in group['raw_files']:
                src = Path(entry['source']).resolve(strict=True)
                src.relative_to(root)
                relative = Path(entry['relative'])
                if relative.is_absolute() or '..' in relative.parts:
                    raise ValueError('Unsafe supplementary path')
                task = src.parent.parent if src.parent.name == 'out' else src.parent
                if (task / '.exitcode').read_text().strip() != '0':
                    raise ValueError('Supplementary task did not finish successfully')
                dst = output / role / relative
                if str(dst) in seen or not src.is_file() or src.stat().st_size <= 0:
                    raise ValueError('Invalid supplementary output')
                seen.add(str(dst))
                rows.append((src, dst, src.stat().st_size, digest(src)))
            continue
        payload = Path(group['payload']).resolve(strict=True)
        payload.relative_to(root)
        if (payload.parent / '.exitcode').read_text().strip() != '0':
            raise ValueError('Payload task did not finish successfully')
        data = json.loads(payload.read_text())
        entries = data.get('files')
        if not isinstance(entries, list) or not entries:
            raise ValueError('Empty payload file list')
        names = [safe_name(e['fileName']) for e in entries]
        if len(set(names)) != len(names):
            raise ValueError('Duplicate payload filename')
        for name in names:
            if name.endswith('.cram') and name + '.crai' not in names:
                raise ValueError('Missing CRAM index')
            if name.endswith('.vcf.gz') and name + '.tbi' not in names:
                raise ValueError('Missing VCF index')
        if role.endswith('pGenDnaAln') and not any(n.endswith('.cram') for n in names):
            raise ValueError('Missing aligned CRAM')
        if role.endswith('pGenDnaSeqQc'):
            for suffix in ('.ubam_qc_metrics.tgz', '.duplicates_metrics.tgz', '.oxog_metrics.tgz', '.qc_metrics.tgz'):
                if not any(n.endswith(suffix) for n in names):
                    raise ValueError('Missing QC category: ' + suffix)
        if role.endswith(('pGenVarSnv', 'pGenVarIndel', 'pGenVarCnv', 'pGenVarSv')):
            if not any(n.endswith('.vcf.gz') for n in names):
                raise ValueError('Missing variant VCF')
        if role.endswith(('pGenQc', 'pGenVarSupp')) and not any(n.endswith('.tgz') for n in names):
            raise ValueError('Missing QC/supplement archive')
        items = [(payload, output / role / payload.name, payload.stat().st_size, digest(payload))]
        for entry in entries:
            src = (payload.parent / 'out' / entry['fileName']).resolve(strict=True)
            src.relative_to(root)
            md5 = entry.get('fileMd5sum', '')
            if not re.fullmatch('[0-9a-fA-F]{32}', md5):
                raise ValueError('Missing/invalid payload MD5')
            items.append((src, output / role / 'out' / entry['fileName'], int(entry['fileSize']), md5.lower()))
        for src, dst, size, md5 in items:
            if not src.is_file() or size <= 0 or src.stat().st_size != size:
                raise ValueError('Missing, empty or size-mismatched source: ' + str(src))
            if str(dst) in seen:
                raise ValueError('Destination collision')
            seen.add(str(dst))
            rows.append((src, dst, size, md5))
    return rows


def publish_file(src, dst, size, md5, force_copy=False):
    src = src.resolve(strict=True)  # Never hard-link a symlink itself.
    initial = stamp(src)
    if initial[2] != size:
        raise RuntimeError('Source size mismatch')
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_symlink():
        raise ValueError('Refusing destination symlink: ' + str(dst))
    old = stamp(dst) if dst.exists() else None
    same = old is not None and os.path.samefile(str(src), str(dst))
    if old and old[2] == size and digest(dst) == md5:
        if stamp(src) != initial or stamp(dst) != old:
            raise RuntimeError('File changed during verification')
        return {'method': 'existing-inode' if same else 'existing-md5', 'bytes': size}
    temp = dst.with_name('.' + dst.name + '.' + uuid.uuid4().hex + '.publish')
    method = 'hardlink'
    try:
        try:
            if force_copy:
                raise OSError(errno.EXDEV, 'Copy requested')
            os.link(str(src), str(temp))
            with temp.open('rb') as linked:
                os.fsync(linked.fileno())
        except OSError as exc:
            if exc.errno not in (errno.EXDEV, errno.EPERM, errno.EACCES, errno.EOPNOTSUPP, errno.ENOSYS, errno.EMLINK):
                raise
            method = 'verified-copy'
            h = hashlib.md5()
            with src.open('rb') as incoming, temp.open('xb') as outgoing:
                for chunk in iter(lambda: incoming.read(8 * 1024 * 1024), b''):
                    outgoing.write(chunk)
                    h.update(chunk)
                outgoing.flush()
                os.fsync(outgoing.fileno())
            if h.hexdigest() != md5 or temp.stat().st_size != size or digest(temp) != md5:
                raise RuntimeError('Copy checksum/size mismatch')
        # Some shared filesystems expose coarse timestamps: hash the linked file
        # after link creation, never rely on size/mtime alone for source stability.
        if method == 'hardlink' and digest(temp) != md5:
            raise RuntimeError('Linked source checksum mismatch')
        # Creating a hard link changes ctime; compare content identity/mtime/size.
        if stamp(src)[:4] != initial[:4] or temp.stat().st_size != size:
            raise RuntimeError('Source changed during publication')
        if (stamp(dst) if dst.exists() else None) != old or dst.is_symlink():
            raise RuntimeError('Destination changed during publication')
        backup = None
        if old:
            backup = dst.with_name('.' + dst.name + '.' + uuid.uuid4().hex + '.previous')
            os.link(str(dst), str(backup))
        os.replace(str(temp), str(dst))  # temp is always on destination filesystem.
        sync_dir(dst.parent)
        if stamp(src)[:4] != initial[:4] or dst.stat().st_size != size:
            raise RuntimeError('Publication verification failed')
        if method == 'hardlink' and not os.path.samefile(str(src), str(dst)):
            raise RuntimeError('Hard-link identity mismatch')
        return {'method': method, 'bytes': size, 'md5': md5, 'backup': str(backup) if backup else None}
    finally:
        if temp.exists():
            temp.unlink()


def publish(manifest, force_copy=False):
    output = Path(manifest['output']).resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = {'version': 1, 'state': 'CHECKING', 'files': []}
    with (output / '.publication.lock').open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        marker = output / 'publication.verified.json'
        if marker.exists():
            marker.unlink()
            sync_dir(output)
        save(marker, report)  # Invalidate an earlier successful run before any change.
        try:
            rows = plan(manifest)  # Validate every required source before publishing.
            for src, dst, size, md5 in rows:
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.parent.resolve().relative_to(output)
                result = publish_file(src, dst, size, md5, force_copy)
                report['files'].append(dict(source=str(src), destination=str(dst), **result))
            report['state'] = 'VERIFIED'
            save(marker, report)
        except Exception as exc:
            report.update(state='FAILED', error=str(exc))
            save(marker, report)
            raise
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest-base64', required=True)
    parser.add_argument('--force-copy', action='store_true')
    args = parser.parse_args()
    manifest = json.loads(base64.b64decode(args.manifest_base64))
    result = publish(manifest, args.force_copy)
    save(Path('publication.verified.json'), result)
    print('VERIFIED: {} files'.format(len(result['files'])))
