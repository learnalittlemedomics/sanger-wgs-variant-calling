nextflow.enable.dsl = 2

// Runs on the Nextflow controller. On clusters, submit that controller to a
// compute allocation. No container mounts or cross-device renames are needed.
process verifiedPublish {
    executor 'local'
    container null
    cpus 1
    memory '1 GB'
    cache false
    errorStrategy 'terminate'
    publishDir enabled: false

    input:
    val groups
    val expected
    val destination

    output:
    path 'publication.verified.json', emit: receipt

    script:
    def target = java.nio.file.Paths.get(destination.toString())
    if (!target.isAbsolute()) target = workflow.launchDir.resolve(target)
    def manifest = [groups: groups, expected: expected, output: target.toString(), work_root: workflow.workDir.toString()]
    def encoded = groovy.json.JsonOutput.toJson(manifest).getBytes('UTF-8').encodeBase64().toString()
    def helper = "${moduleDir}/../scripts/verified_publish.py".replace("'", "'\"'\"'")
    """
    python3 '${helper}' --manifest-base64 '${encoded}'
    """
}
