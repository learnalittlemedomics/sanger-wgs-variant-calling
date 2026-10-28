nextflow.enable.dsl=2
process SubmitRetryProbe {
    cpus 1
    memory '64 MB'
    errorStrategy 'terminate'
    output:
    path 'attempt.txt'
    script:
    """
    echo '${task.attempt}' > attempt.txt
    """
}
workflow { SubmitRetryProbe() }
