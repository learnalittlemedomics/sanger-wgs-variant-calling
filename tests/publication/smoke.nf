nextflow.enable.dsl = 2
include { verifiedPublish } from '../../modules/verified-publish'
params.outdir = "${launchDir}/smoke-output"
params.broken = false
process payload {
    output:
    path 'test.payload.json', emit: payload
    path 'out/*', emit: files
    script:
    """
    mkdir out
    python3 - <<'PY'
    import json,hashlib
    from pathlib import Path
    files=[]
    for n in ['call.vcf.gz','call.vcf.gz.tbi']:
        p=Path('out')/n;p.write_bytes(n.encode()*100)
        files.append(dict(fileName=n,fileSize=p.stat().st_size,fileMd5sum=hashlib.md5(p.read_bytes()).hexdigest()))
    Path('test.payload.json').write_text(json.dumps(dict(files=files)))
    if '${params.broken}' == 'true': Path('out/call.vcf.gz.tbi').unlink()
    PY
    """
}
process afterPublication {
    input:
    path receipt
    output:
    path 'after.txt'
    script:
    """
    python3 -c 'import json; assert json.load(open("${receipt}"))["state"] == "VERIFIED"'
    touch after.txt
    """
}
workflow {
    payload()
    groups=payload.out.payload.map { [role:'M2_pGenVarSnv', payload:it.toString()] }
    verifiedPublish(groups.collect(flat:false), ['M2_pGenVarSnv'], params.outdir)
    afterPublication(verifiedPublish.out.receipt)
}
