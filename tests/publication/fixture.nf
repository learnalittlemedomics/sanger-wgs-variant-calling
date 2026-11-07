nextflow.enable.dsl=2
process payload {
    input:
    val role
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
    role='${role}'
    names=['call.vcf.gz','call.vcf.gz.tbi']
    if role.endswith('pGenDnaAln'): names=['aln.cram','aln.cram.crai']
    if role.endswith('pGenDnaSeqQc'): names=['x'+s for s in ['.ubam_qc_metrics.tgz','.duplicates_metrics.tgz','.oxog_metrics.tgz','.qc_metrics.tgz']]
    if role.endswith(('pGenQc','pGenVarSupp')): names=['metrics.tgz']
    for n in names:
        p=Path('out')/n;p.write_bytes(n.encode()*100)
        files.append(dict(fileName=n,fileSize=p.stat().st_size,fileMd5sum=hashlib.md5(p.read_bytes()).hexdigest()))
    Path('test.payload.json').write_text(json.dumps(dict(files=files)))
    
    PY
    """
}
