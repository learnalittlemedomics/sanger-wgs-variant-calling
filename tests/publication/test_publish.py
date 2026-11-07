import sys, tempfile, json, hashlib, errno, os
from pathlib import Path
from unittest.mock import patch
here=Path(__file__).resolve().parent
scripts=here/'scripts' if (here/'scripts').is_dir() else here.parent.parent/'scripts'
sys.path.insert(0,str(scripts))
import verified_publish as v

def fixture(root):
    task=root/'work'/'ab'/'cdef';task.mkdir(parents=True)
    (task/'.exitcode').write_text('0');(task/'out').mkdir()
    source=root/'work'/'source';source.mkdir()
    files=[]
    for name in ['call.vcf.gz','call.vcf.gz.tbi']:
        p=source/name;p.write_bytes((name*100).encode())
        (task/'out'/name).symlink_to(p)
        files.append(dict(fileName=name,fileSize=p.stat().st_size,fileMd5sum=v.digest(p)))
    payload=task/'test.payload.json';payload.write_text(json.dumps(dict(files=files)))
    return dict(work_root=str(root/'work'),output=str(root/'output'),expected=['M2_pGenVarSnv'],groups=[dict(role='M2_pGenVarSnv',payload=str(payload))])

def test():
    with tempfile.TemporaryDirectory(dir='.') as folder:
        root=Path(folder).resolve();m=fixture(root)
        r=v.publish(m);assert r['state']=='VERIFIED' and len(r['files'])==3
        dst=root/'output/M2_pGenVarSnv/out/call.vcf.gz.tbi';src=root/'work/source/call.vcf.gz.tbi'
        assert not dst.is_symlink() and os.path.samefile(str(dst),str(src))
        # Replace (not truncate) a hard link to keep the source immutable.
        dst.unlink();dst.write_bytes(b'bad')
        v.publish(m);assert os.path.samefile(str(dst),str(src))
        assert list(dst.parent.glob('.*.previous'))
        # Equal-size corruption must also be repaired.
        dst.unlink();dst.write_bytes(b'x'*src.stat().st_size);v.publish(m)
        # Missing required output invalidates a previous success.
        src.rename(src.with_suffix('.saved'))
        try:v.publish(m);raise AssertionError('Missing index accepted')
        except FileNotFoundError:pass
        assert json.loads((root/'output/publication.verified.json').read_text())['state']=='FAILED'
    with tempfile.TemporaryDirectory(dir='.') as folder:
        root=Path(folder).resolve();m=fixture(root)
        v.publish(m,force_copy=True)
        dst=root/'output/M2_pGenVarSnv/out/call.vcf.gz';src=root/'work/source/call.vcf.gz'
        assert not os.path.samefile(str(dst),str(src)) and v.digest(dst)==v.digest(src)
    with tempfile.TemporaryDirectory(dir='.') as folder:
        root=Path(folder).resolve();m=fixture(root)
        real_link=os.link
        def exdev(src,dst,*a,**kw):
            if str(dst).endswith('.publish'):raise OSError(errno.EXDEV,'test')
            return real_link(src,dst,*a,**kw)
        with patch.object(v.os,'link',side_effect=exdev):r=v.publish(m)
        assert all(f['method']=='verified-copy' for f in r['files'])
    with tempfile.TemporaryDirectory(dir='.') as folder:
        root=Path(folder).resolve();m=fixture(root)
        v.publish(m,force_copy=True)
        dst=root/'output/M2_pGenVarSnv/out/call.vcf.gz.tbi'
        dst.write_bytes(b'bad')
        with patch.object(v.os,'link',side_effect=OSError(errno.EOPNOTSUPP,'Unsupported')):
            v.publish(m)
        assert dst.read_bytes()==(root/'work/source/call.vcf.gz.tbi').read_bytes()
        assert any(p.read_bytes()==b'bad' for p in dst.parent.glob('.*.previous'))
    with tempfile.TemporaryDirectory(dir='.') as folder:
        root=Path(folder).resolve();m=fixture(root)
        # Failed copy/checksum cannot expose a final file or leave VERIFIED.
        payload=Path(m['groups'][0]['payload']);d=json.loads(payload.read_text());d['files'][0]['fileMd5sum']='0'*32;payload.write_text(json.dumps(d))
        try:v.publish(m,force_copy=True);raise AssertionError('Bad checksum accepted')
        except RuntimeError:pass
        assert not (root/'output/M2_pGenVarSnv/out/call.vcf.gz').exists()
        assert json.loads((root/'output/publication.verified.json').read_text())['state']=='FAILED'
    with tempfile.TemporaryDirectory(dir='.') as folder:
        root=Path(folder).resolve();m=fixture(root)
        src=root/'work/source/call.vcf.gz'
        real_link=os.link
        def changing(source,destination,*a,**kw):
            result=real_link(source,destination,*a,**kw)
            if str(source)==str(src):src.write_bytes(b'x'*src.stat().st_size)
            return result
        with patch.object(v.os,'link',side_effect=changing):
            try:v.publish(m);raise AssertionError('Changing source accepted')
            except RuntimeError:pass
        assert not (root/'output/M2_pGenVarSnv/out/call.vcf.gz').exists()
    with tempfile.TemporaryDirectory(dir='.') as folder:
        root=Path(folder).resolve();m=fixture(root)
        # ENOSPC must not be treated as an unsupported-link fallback.
        with patch.object(v.os,'link',side_effect=OSError(errno.ENOSPC,'Full')):
            try:v.publish(m);raise AssertionError('Full disk accepted')
            except OSError as exc:assert exc.errno==errno.ENOSPC
        assert json.loads((root/'output/publication.verified.json').read_text())['state']=='FAILED'
    # True cross-device publication on the compute node: ParaStor -> ext4 /tmp.
    with tempfile.TemporaryDirectory(dir='.') as folder, tempfile.TemporaryDirectory(dir='/tmp') as other:
        root=Path(folder).resolve();m=fixture(root);m['output']=other
        if root.stat().st_dev!=Path(other).stat().st_dev:
            r=v.publish(m);assert all(f['method']=='verified-copy' for f in r['files'])
            print('REAL_CROSS_DEVICE_PASSED',root.stat().st_dev,Path(other).stat().st_dev)
        else:
            assert not os.environ.get('REQUIRE_CROSS_DEVICE'), 'Need distinct mounts'
            print('REAL_CROSS_DEVICE_SKIPPED: /tmp and current directory share a device')
    print('PYTHON_TESTS_PASSED')

if __name__=='__main__':test()
