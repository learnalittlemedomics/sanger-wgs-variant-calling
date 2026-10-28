nextflow.enable.dsl = 2
params.cpus = 24
params.mem = 128
params.container_version = ''
params.exclude = 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV'
params.species = 'human'
params.assembly = 'GRCh38'
params.skipqc = false
params.skipannot = true
params.pindelcpu = 8
params.cavereads = 800000
params.purity = null
params.ploidy = null
params.steps = [:]
params.split_sif = ""
params.resources = [:]
params.scatter = [:]
params.dependency_optimize = true


// Only files actually executed by this stage enter its Nextflow input hash.
def stageScripts(root, stage) {
    def names
    if (stage in ['CaVEManFlagPlan','CaVEManFlagGroup','CaVEMan_flag']) names = ['run-flag-scatter.sh','flag-scatter.py','command-CaVEMan_flag.sh']
    else if (stage == 'BRASSGather') names = ['run-brass-gather.sh','brass-gather.py','command-BRASS.sh']
    else if (stage.startsWith('BRASS') && !(stage in ['BRASS','BRASS_input','BRASS_cover'])) names = ['run-brass-scatter.sh','brass-scatter.py','command-BRASS.sh']
    else if (stage == 'prepareEarly') names = ['prepare-early.pl']
    else if (stage == 'cacheEarly') names = ['cache-early.sh']
    else if (stage == 'PindelInputEarly') names = ['run-pindel-early16.sh','command-cgpPindel.sh','pindel-input16.py','pindel-layout.pl','pindel_input_gen16.pl']
    else if (stage == 'PindelInput') names = ['run-pindel-input16.sh','command-cgpPindel.sh','pindel-input16.py','pindel-layout.pl','pindel_input_gen16.pl']
    else if (stage == 'PindelFinalize') names = ['run-pindel-finalize.sh','command-cgpPindel.sh','pindel-finalize.py','pindel-layout.pl']
    else if (stage.startsWith('Pindel')) names = ['run-pindel-scatter.sh','command-cgpPindel.sh','pindel-scatter.py','pindel-layout.pl']
    else if (stage == 'prepare') names = ['prepare.pl']
    else if (stage == 'assemble') names = ['assemble.sh']
    else if (stage in ['CaVEMan_Mstep','CaVEMan_ModelMerge','CaVEMan_Estep','CaVEMan_Gather'])
        names = ['caveman-scatter.py']
    else {
        names = ['run-stage.sh']
        if (stage != 'cache_POP')
            names += ['command-' + (stage == 'CaVEMan_plan' ? 'CaVEMan' : stage) + '.sh']
        if (stage == 'CaVEMan_plan') names += ['caveman-scatter.py']
    }
    if (stage == 'CaVEManFlagGroup') names += ['CavemanPostProcessor.pm']
    if (stage == 'BRASSPlan') names += ['get_abs_bkpts_from_clipped_reads.pl']
    return names.collect { root.resolve(it) }
}

def bindInput16(task, engine) {
    def src = "${task.workDir}/split_scripts/pindel_input_gen16.pl"
    def dst = '/opt/wtsi-cgp/bin/pindel_input_gen.pl'
    def flag = engine in ['singularity','apptainer'] ? '--bind' : '-v'
    return bindWork(task, engine) + " ${flag} ${src}:${dst}:ro"
}

def bindWork(task, engine) {
    if (engine in ['singularity', 'apptainer']) return "--bind ${task.workDir}:/sanger"
    if (engine in ['docker', 'podman']) return "-v ${task.workDir}:/sanger"
    throw new IllegalArgumentException('Split Sanger requires Singularity/Apptainer or Docker/Podman for stable /sanger paths')
}

// Each allocation runs a bounded pool of single-threaded CaVEMan lanes.
def caveCpus(task, kind) {
    def defaults = kind == 'group' ? [30,30,30,30] : [1,2,3,5]
    def steps = params.scatter[kind + '_cpu_steps'] ?: defaults
    if (!(steps instanceof List) || !steps || steps.any { !(it instanceof Number) || it.intValue() != it || it < 1 || it > 30 })
        throw new IllegalArgumentException('scatter CPU steps must be integers in 1..30')
    return steps[Math.min((task.attempt as int)-1, steps.size()-1)] as int
}
def caveWorkers(task) {
    def steps = params.scatter.worker_steps ?: [12,8,6,4]
    def packing = (params.scatter.workers_per_job ?: 12) as int
    if (!(steps instanceof List) || !steps || packing < 1 || packing > 30 ||
        steps.any { !(it instanceof Number) || it.intValue() != it || it < 1 || it > packing })
        throw new IllegalArgumentException('worker_steps must be integers in 1..workers_per_job (<=30)')
    def n = steps[Math.min((task.attempt as int)-1, steps.size()-1)] as int
    if (n > caveCpus(task, 'group')) throw new IllegalArgumentException('workers exceed allocated CPUs')
    return n
}
def caveMemory(task, kind) {
    def gb = (params.scatter.memory_gb_per_cpu ?: 3.6) as double
    if (gb <= 0 || gb > 3.6) throw new IllegalArgumentException('memory_gb_per_cpu must be in (0,3.6] for this cluster')
    // Round down in MiB so the request never exceeds the per-core allocation.
    return "${Math.floor(caveCpus(task, kind) * gb * 1024) as long} MB"
}
def caveRetry(task) {
    def code = task.exitStatus
    // 2147483647: Nextflow reports unknown external termination / missing .exitcode.
    if (code in [137, 140, 143, 247, 2147483647] || code == null) return 'retry'
    for (name in ['.command.err', '.command.log']) {
        def logFile = task.workDir.resolve(name)
        if (logFile.exists()) {
            def text = logFile.text
            if (text =~ /(?i)(out.of.memory|oom[_ -]?kill|cannot allocate memory|MemoryError|SCATTER_RESOURCE_FAILURE)/) return 'retry'
        }
    }
    return 'terminate'
}
process prepare {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_prepare'
    input:
    path reference
    path annot
    path snv_indel
    path cnv_sv
    path qcset
    path tumour
    path tidx
    path tumour_bas
    path normal
    path nidx
    path normal_bas
    path scripts, stageAs: 'split_scripts/*'
    output:
    tuple path('reference_files'), path('run.params'), path('flag.vcf.config.WGS.ini'), path(tumour), path(ti), path(tumour_bas), path(normal), path(ni), path(normal_bas), emit: context
    script:
    def tlist = tidx instanceof List ? tidx : [tidx]
    def nlist = nidx instanceof List ? nidx : [nidx]
    ti = tlist.find { it.name.endsWith(tumour.name.endsWith('.cram') ? '.crai' : '.bai') } ?: tlist[0]
    ni = nlist.find { it.name.endsWith(normal.name.endsWith('.cram') ? '.crai' : '.bai') } ?: nlist[0]
    def pu = params.purity != null ? "-pu ${params.purity}" : ''
    def pi = params.ploidy != null ? "-pi ${params.ploidy}" : ''
    """
    cd /sanger
    perl split_scripts/prepare.pl -cores ${task.cpus} \
      -reference ${reference} -annot ${annot} -snv_indel ${snv_indel} \
      -cnv_sv ${cnv_sv} -qcset ${qcset} \
      -tumour ${tumour} -tidx ${ti} -normal ${normal} -nidx ${ni} \
      -exclude '${params.exclude}' -species '${params.species}' -assembly '${params.assembly}' \
      ${params.skipqc ? '-skipqc' : ''} ${params.skipannot ? '-skipannot' : ''} \
      -pindelcpu ${params.pindelcpu} -cavereads ${params.cavereads} ${pu} ${pi} -outdir /sanger
    """
}
process cache_POP {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_cache_POP'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    output:
    path 'cache_POP.stage.tar.gz', emit: archive
    path 'ref_cache', emit: cache
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh cache_POP ${task.cpus}
    """
}
process CaVEMan_setup {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_CaVEMan_setup'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'CaVEMan_setup.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh CaVEMan_setup ${task.cpus} ${dependencies}
    """
}
process geno {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_geno'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'geno.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh geno ${task.cpus} ${dependencies}
    """
}
process verify_WT {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_verify_WT'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'verify_WT.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh verify_WT ${task.cpus} ${dependencies}
    """
}
process CaVEMan_split {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_CaVEMan_split'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'CaVEMan_split.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh CaVEMan_split ${task.cpus} ${dependencies}
    """
}
process ascat {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_ascat'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'ascat.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh ascat ${task.cpus} ${dependencies}
    """
}
process BRASS_input {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_BRASS_input'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS_input.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh BRASS_input ${task.cpus} ${dependencies}
    """
}
process BRASS_cover {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_BRASS_cover'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS_cover.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh BRASS_cover ${task.cpus} ${dependencies}
    """
}
process cgpPindel {
    label 'sanger_split'
    cpus 30
    memory '108 GB'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_cgpPindel'
    maxRetries 2
    errorStrategy { task.exitStatus in [137, 140, 143, 247] ? 'retry' : 'terminate' }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'cgpPindel.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    export PINDEL_WORKERS='${task.attempt == 1 ? 16 : task.attempt == 2 ? 12 : 8}'
    bash split_scripts/run-stage.sh cgpPindel ${task.cpus} ${dependencies}
    """
}
process PindelInput {
    label 'sanger_split'
    label 'sanger_pindel_scatter'
    maxRetries 3
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindInput16(task, workflow.containerEngine) }
    cpus { [20,24,30,30][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    val index
    output:
    tuple val(index), path('pindel_input_*'), emit: parts
    script:
    """
    cd /sanger
    bash split_scripts/run-pindel-input16.sh input ${task.cpus} --index ${index}
    """
    stub:
    """
    mkdir -p pindel_input_${index}/data; echo ok > pindel_input_${index}/data/mock
    """
}
process PindelPlan {
    label 'sanger_split'
    label 'sanger_pindel_scatter'
    maxRetries 3
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    cpus { [1,2,3,5][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path parts, stageAs: 'pindel_inputs/*'
    output:
    path 'pindel-plan.json', emit: plan
    script:
    """
    cd /sanger
    bash split_scripts/run-pindel-scatter.sh plan ${task.cpus} --inputs ${parts}
    """
    stub:
    """
    echo '{"jobs":2,"chromosomes":["chr1","chr2"]}' > pindel-plan.json
    """
}
process PindelChromosome {
    label 'sanger_split'
    label 'sanger_pindel_scatter'
    maxRetries 3
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    cpus { [8,12,20,30][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    tag "${chromosome}"
    maxForks 24
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path parts, stageAs: 'pindel_inputs/*'
    path plan, stageAs: 'pindel-plan.json'
    tuple val(index), val(chromosome)
    output:
    tuple val(index), path('pindel_part_*'), emit: parts
    script:
    """
    cd /sanger
    bash split_scripts/run-pindel-scatter.sh chromosome ${task.cpus} --index ${index} --inputs ${parts} --plan ${plan}
    """
    stub:
    """
    if [ "${index}" = 1 ] && [ "${task.attempt}" = 1 ]; then exit 137; fi
    mkdir -p pindel_part_${index}/vcf; echo ok > pindel_part_${index}/vcf/mock
    """
}
process PindelGather {
    label 'sanger_split'
    label 'sanger_pindel_scatter'
    maxRetries 3
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    cpus { [8,12,20,30][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path plan, stageAs: 'pindel-plan.json'
    path parts, stageAs: 'pindel_parts/*'
    output:
    path 'cgpPindel.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-pindel-scatter.sh gather ${task.cpus} --plan ${plan} --parts ${parts}
    """
    stub:
    """
    touch cgpPindel.stage.tar.gz
    """
}

process CaVEMan_plan {
    label 'sanger_split'
    label 'sanger_caveman_scatter'
    cpus { caveCpus(task, 'small') }
    memory { caveMemory(task, 'small') }
    maxRetries { (params.scatter.max_retries == null ? 3 : params.scatter.max_retries) as int }
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'scatter_context', emit: context
    path 'scatter-plan.json', emit: plan
    script:
    """
    cd /sanger
    export CAVE_WORKERS_PER_JOB='${params.scatter.workers_per_job ?: 12}'
    export CAVE_TARGET_HOURS='${params.scatter.target_hours ?: 1.5}'
    export CAVE_MSTEP_MINUTES='${params.scatter.mstep_minutes_per_shard ?: 3.328}'
    export CAVE_ESTEP_MINUTES='${params.scatter.estep_minutes_per_shard ?: 15.83}'
    bash split_scripts/run-stage.sh CaVEMan_scatter_prepare 1 ${dependencies}
    """
    stub:
    """
    mkdir -p scatter_context
    echo '{"mstep_groups":2,"estep_groups":3}' > scatter-plan.json
    cp scatter-plan.json scatter_context/plan.json
    """
}

process CaVEMan_Mstep {
    label 'sanger_split'
    label 'sanger_caveman_scatter'
    cpus { caveCpus(task, 'group') }
    memory { caveMemory(task, 'group') }
    maxRetries { (params.scatter.max_retries == null ? 3 : params.scatter.max_retries) as int }
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    tag "group_${index}"
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    tuple val(index), path(context)
    output:
    tuple val(index), path('result'), emit: parts
    script:
    """
    cd /sanger
    python3 split_scripts/caveman-scatter.py run --mode mstep --index ${index} --workers ${caveWorkers(task)}
    """
    stub:
    """
    if [ "${index}" = 1 ] && [ "${task.attempt}" -le 3 ]; then exit 137; fi
    mkdir -p result
    echo 'attempt=${task.attempt} cpus=${task.cpus} memory=${task.memory} workers=${caveWorkers(task)}' > result/resources.txt
    """
}

process CaVEMan_ModelMerge {
    label 'sanger_split'
    label 'sanger_caveman_scatter'
    cpus { caveCpus(task, 'small') }
    memory { caveMemory(task, 'small') }
    maxRetries { (params.scatter.max_retries == null ? 3 : params.scatter.max_retries) as int }
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path context, stageAs: 'scatter_context'
    path parts, stageAs: 'parts/part????/*'
    output:
    path 'result', emit: model
    script:
    def quoted = parts.collect { "'${it}'" }.join(' ')
    """
    cd /sanger
    python3 split_scripts/caveman-scatter.py run --mode merge --parts ${quoted}
    """
    stub:
    """
    mkdir -p result
    echo 'attempt=${task.attempt} cpus=${task.cpus} memory=${task.memory}' > result/resources.txt
    """
}

process CaVEMan_Estep {
    label 'sanger_split'
    label 'sanger_caveman_scatter'
    cpus { caveCpus(task, 'group') }
    memory { caveMemory(task, 'group') }
    maxRetries { (params.scatter.max_retries == null ? 3 : params.scatter.max_retries) as int }
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    tag "group_${index}"
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    tuple val(index), path(context)
    path model, stageAs: 'model'
    output:
    tuple val(index), path('result'), emit: parts
    script:
    """
    cd /sanger
    python3 split_scripts/caveman-scatter.py run --mode estep --index ${index} --workers ${caveWorkers(task)}
    """
    stub:
    """
    mkdir -p result
    echo 'attempt=${task.attempt} cpus=${task.cpus} memory=${task.memory} workers=${caveWorkers(task)}' > result/resources.txt
    """
}

process CaVEMan_Gather {
    label 'sanger_split'
    label 'sanger_caveman_scatter'
    cpus { caveCpus(task, 'small') }
    memory { caveMemory(task, 'small') }
    maxRetries { (params.scatter.max_retries == null ? 3 : params.scatter.max_retries) as int }
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path context, stageAs: 'scatter_context'
    path parts, stageAs: 'parts/part????/*'
    path model, stageAs: 'model'
    output:
    path 'CaVEMan.stage.tar.gz', emit: archive
    path 'result/report.json', emit: report
    script:
    def quoted = parts.collect { "'${it}'" }.join(' ')
    """
    cd /sanger
    python3 split_scripts/caveman-scatter.py run --mode gather --parts ${quoted}
    """
    stub:
    """
    mkdir -p result
    echo PASS > result/report.json
    tar -czf CaVEMan.stage.tar.gz result
    """
}

process BRASSInput30 {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    cpus 30
    memory '110592 MB'
    tag "input_${role}"
    input:
    val role
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    tuple val(role), path("BRASS_input_${role}.stage.tar"), emit: parts
    script:
    """
    cd /sanger
    export BRASS_GROUP='${role}'
    bash split_scripts/run-brass-scatter.sh input ${task.cpus} ${dependencies}
    """
    stub:
    """
    touch BRASS_input_${role}.stage.tar
    """
}

process BRASSCoverage {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    cpus 30
    memory '110592 MB'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS_cover.stage.tar', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-brass-scatter.sh cover ${task.cpus} ${dependencies}
    """
    stub:
    """
    touch BRASS_cover.stage.tar
    """
}

process BRASSGroup {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    cpus 30
    memory '110592 MB'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS_group.stage.tar', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-brass-scatter.sh group ${task.cpus} ${dependencies}
    """
    stub:
    """
    touch BRASS_group.stage.tar
    """
}

process BRASSInsertSize {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    cpus 4
    memory '14745 MB'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS_isize.stage.tar', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-brass-scatter.sh isize ${task.cpus} ${dependencies}
    """
    stub:
    """
    touch BRASS_isize.stage.tar
    """
}

process BRASSNormCN {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    cpus { [2,4,8,16][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS_normcn.stage.tar', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-brass-scatter.sh normcn ${task.cpus} ${dependencies}
    """
    stub:
    """
    touch BRASS_normcn.stage.tar
    """
}

process BRASSPlan {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions {
        def src = "${task.workDir}/split_scripts/get_abs_bkpts_from_clipped_reads.pl"
        def dst = '/opt/wtsi-cgp/bin/get_abs_bkpts_from_clipped_reads.pl'
        def flag = workflow.containerEngine in ['singularity','apptainer'] ? '--bind' : '-v'
        return bindWork(task, workflow.containerEngine) + " ${flag} ${src}:${dst}:ro"
    }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    cpus 30
    memory '110592 MB'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS_plan.stage.tar', emit: archive
    path 'brass-seed.tar', emit: seed
    path 'brass-plan.json', emit: plan
    script:
    """
    cd /sanger
    bash split_scripts/run-brass-scatter.sh plan ${task.cpus} ${dependencies}
    """
    stub:
    """
    touch BRASS_plan.stage.tar
    touch brass-seed.tar
    echo '{"splits":7,"groups":[[1,2,3,4,5,6],[7]]}' > brass-plan.json
    """
}

process BRASSGather {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    cpus { [2,4,8,16][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-brass-gather.sh gather ${task.cpus} ${dependencies}
    """
    stub:
    """
    touch BRASS.stage.tar.gz
    """
}

process BRASSAssemble {
    label 'sanger_split'
    label 'sanger_brass_scatter'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    errorStrategy { caveRetry(task) }
    maxRetries 3
    tag "group_${group}"
    cpus { [24,30,30,30][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    maxForks 32
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    tuple val(group), path(seed)
    output:
    tuple val(group), path("brass-part-${group}.tar"), emit: parts
    script:
    def workers = [6,6,4,2][Math.min(task.attempt-1,3)]
    """
    cd /sanger
    export BRASS_GROUP='${group}' BRASS_WORKERS='${workers}'
    bash split_scripts/run-brass-scatter.sh assemble ${task.cpus} ${seed}
    """
    stub:
    """
    touch brass-part-${group}.tar
    """
}

process BRASS {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_BRASS'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'BRASS.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh BRASS ${task.cpus} ${dependencies}
    """
}
process cgpPindel_annot {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_cgpPindel_annot'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'cgpPindel_annot.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh cgpPindel_annot ${task.cpus} ${dependencies}
    """
}
process CaVEManFlagPlan {
    label 'sanger_split'
    cpus 1
    memory '3686 MB'
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'flag-plan.tar.gz', emit: archive
    path 'flag-groups/*.tar.gz', emit: splits
    path 'flag-options.json', emit: options
    script:
    """
    cd /sanger
    FLAG_MODE=plan bash split_scripts/run-flag-scatter.sh CaVEMan_flag ${task.cpus} ${dependencies}
    """
}
process CaVEManFlagGroup {
    tag "group_$shard_id"
    maxForks ((params.scatter.flag_max_forks ?: 8) as int)
    errorStrategy { task.exitStatus in [137,140,143] ? 'retry' : 'terminate' }
    maxRetries 2
    label 'sanger_split'
    cpus 30
    memory { "${Math.floor(3.6 * 1024 * 30) as long} MB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions {
        def src = "${task.workDir}/split_scripts/CavemanPostProcessor.pm"
        def dst = '/opt/wtsi-cgp/lib/perl5/Sanger/CGP/CavemanPostProcessor.pm'
        def flag = workflow.containerEngine in ['singularity','apptainer'] ? '--bind' : '-v'
        return bindWork(task, workflow.containerEngine) + " ${flag} ${src}:${dst}:ro"
    }
    
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    path flag_options, stageAs: 'flag-options.json'
    tuple val(shard_id), path(split_vcf)
    output:
    path 'flag-part-*.tar.gz', emit: archive
    script:
    """
    cd /sanger
    FLAG_MODE=worker FLAG_WORKERS=${Math.min(task.cpus as int, task.attempt == 1 ? 30 : task.attempt == 2 ? 15 : 10)} FLAG_ID=${shard_id} FLAG_INPUT=${split_vcf} bash split_scripts/run-flag-scatter.sh CaVEMan_flag ${task.cpus} ${dependencies}
    """
}
process CaVEMan_flag {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_CaVEMan_flag'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'CaVEMan_flag.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    FLAG_MODE=gather bash split_scripts/run-flag-scatter.sh CaVEMan_flag ${task.cpus} ${dependencies}
    """
}
process CaVEMan_annot {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_CaVEMan_annot'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'CaVEMan_annot.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh CaVEMan_annot ${task.cpus} ${dependencies}
    """
}
process verify_MT {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_verify_MT'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path dependencies, stageAs: 'deps/*'
    output:
    path 'verify_MT.stage.tar.gz', emit: archive
    script:
    """
    cd /sanger
    bash split_scripts/run-stage.sh verify_MT ${task.cpus} ${dependencies}
    """
}
process assemble {
    label 'sanger_split'
    cpus { params.resources[task.process.tokenize(':')[-1]]?.cpus ?: params.cpus }
    memory { "${params.resources[task.process.tokenize(':')[-1]]?.mem ?: params.mem} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_assemble'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path archives, stageAs: 'archives/*'
    output:
    path 'run.params', emit: run_params
    path 'WGS_*_vs_*.result.tar.gz', emit: result_archive
    path 'WGS_*_vs_*.timings.tar.gz', emit: timings
    path 'WGS_*_vs_*.time', emit: global_time
    script:
    """
    cd /sanger
    bash split_scripts/assemble.sh ${archives}
    """
}
process prepareEarly {
    label 'sanger_split'
    cpus { params.resources['prepare']?.cpus ?: 5 }
    memory { "${params.resources['prepare']?.mem ?: 18} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_prepare'
    input:
    path reference
    path annot
    path snv_indel
    path cnv_sv
    path qcset
    path tumour
    path tidx
    path normal
    path nidx
    path scripts, stageAs: 'split_scripts/*'
    output:
    tuple path('reference_files'), path('run.params'), path('flag.vcf.config.WGS.ini'), path(tumour), path(ti), path(normal), path(ni), emit: context
    script:
    def tlist = tidx instanceof List ? tidx : [tidx]
    def nlist = nidx instanceof List ? nidx : [nidx]
    ti = tlist.find { it.name.endsWith(tumour.name.endsWith('.cram') ? '.crai' : '.bai') } ?: tlist[0]
    ni = nlist.find { it.name.endsWith(normal.name.endsWith('.cram') ? '.crai' : '.bai') } ?: nlist[0]
    def pu = params.purity != null ? "-pu ${params.purity}" : ''
    def pi = params.ploidy != null ? "-pi ${params.ploidy}" : ''
    """
    cd /sanger
    perl split_scripts/prepare-early.pl -cores ${task.cpus} \
      -reference ${reference} -annot ${annot} -snv_indel ${snv_indel} \
      -cnv_sv ${cnv_sv} -qcset ${qcset} \
      -tumour ${tumour} -tidx ${ti} -normal ${normal} -nidx ${ni} \
      -exclude '${params.exclude}' -species '${params.species}' -assembly '${params.assembly}' \
      ${params.skipqc ? '-skipqc' : ''} ${params.skipannot ? '-skipannot' : ''} \
      -pindelcpu ${params.pindelcpu} -cavereads ${params.cavereads} ${pu} ${pi} -outdir /sanger
    """
    stub:
    def ts = tidx instanceof List ? tidx : [tidx]
    def ns = nidx instanceof List ? nidx : [nidx]
    ti = ts[0]
    ni = ns[0]
    """
    mkdir -p reference_files
    touch reference_files/genome.fa run.params flag.vcf.config.WGS.ini
    """
}
process cacheEarly {
    label 'sanger_split'
    cpus { params.resources['cache_POP']?.cpus ?: 5 }
    memory { "${params.resources['cache_POP']?.mem ?: 18} GB" }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    label 'sanger_cache_POP'
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(normal), path(nidx)
    path scripts, stageAs: 'split_scripts/*'
    output:
    path 'cache_POP.stage.tar.gz', emit: archive
    path 'ref_cache', emit: cache
    script:
    """
    cd /sanger
    bash split_scripts/cache-early.sh
    """
    stub:
    """
    mkdir -p ref_cache
    touch ref_cache/mock cache_POP.stage.tar.gz
    """
}
process PindelInputEarly {
    label 'sanger_split'
    label 'sanger_pindel_scatter'
    maxRetries 3
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindInput16(task, workflow.containerEngine) }
    cpus { [20,24,30,30][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(normal), path(nidx)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    tuple val(index), path(own_bas)
    output:
    tuple val(index), path('pindel_input_*'), emit: parts
    script:
    """
    cd /sanger
    bash split_scripts/run-pindel-early16.sh input ${task.cpus} --index ${index}
    """
    stub:
    """
    mkdir -p pindel_input_${index}/data; echo ok > pindel_input_${index}/data/mock
    """
}
process PindelFinalize {
    label 'sanger_split'
    label 'sanger_pindel_scatter'
    maxRetries 3
    maxErrors -1
    errorStrategy { caveRetry(task) }
    container { params.split_sif ?: "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: '2.1.0-9'}" }
    containerOptions { bindWork(task, workflow.containerEngine) }
    cpus { [8,12,20,30][Math.min(task.attempt-1,3)] }
    memory { "${Math.floor(task.cpus * 3.6 * 1024) as long} MB" }
    input:
    tuple path(reference_files), path(run_params), path(flag_ini), path(tumour), path(tidx), path(tumour_bas), path(normal), path(nidx), path(normal_bas)
    path scripts, stageAs: 'split_scripts/*'
    path cache, stageAs: 'ref_cache'
    path plan, stageAs: 'pindel-plan.json'
    path parts, stageAs: 'pindel_parts/*'
    output:
    path 'pindel-germline.stage.tar.gz', emit: germline
    tuple path('WGS_*'), path('timings'), emit: result
    script:
    """
    cd /sanger
    bash split_scripts/run-pindel-finalize.sh gather ${task.cpus} --plan ${plan} --parts ${parts}
    """
    stub:
    """
    mkdir -p WGS_mock/pindel timings
    touch pindel-germline.stage.tar.gz
    """
}

process PindelPack {
    label 'sanger_split'
    cpus 2
    memory '7 GB'
    input:
    tuple path(pair), path(timings)
    output:
    path 'cgpPindel.stage.tar.gz', emit: archive
    script:
    """
    pair_dir=\$(readlink -f '${pair}')
    timings_dir=\$(readlink -f '${timings}')
    tar -czf cgpPindel.stage.tar.gz -C "\$(dirname "\$pair_dir")" "\$(basename "\$pair_dir")" -C "\$(dirname "\$timings_dir")" "\$(basename "\$timings_dir")"
    """
    stub:
    """
    touch cgpPindel.stage.tar.gz
    """
}

workflow sangerWgsSplit {
    take:
    reference
    annot
    snv_indel
    cnv_sv
    qcset
    tumour
    tidx
    tumour_bas
    normal
    nidx
    normal_bas
    scripts
    main:
    if (params.dependency_optimize) {
        prepareEarly(reference, annot, snv_indel, cnv_sv, qcset, tumour, tidx, normal, nidx, stageScripts(scripts, 'prepareEarly'))
        raw_ctx = prepareEarly.out.context.first()
        ctx = raw_ctx.combine(tumour_bas).combine(normal_bas).map { r, p, f, t, ti, n, ni, tb, nb -> tuple(r,p,f,t,ti,tb,n,ni,nb) }.first()
        cacheEarly(raw_ctx, stageScripts(scripts, 'cacheEarly'))
        cache = cacheEarly.out.cache.first()
        reference_archive = cacheEarly.out.archive
    } else {
        prepare(reference, annot, snv_indel, cnv_sv, qcset, tumour, tidx, tumour_bas, normal, nidx, normal_bas, stageScripts(scripts, 'prepare'))
        ctx = prepare.out.context.first()
        cache_POP(ctx, stageScripts(scripts, 'cache_POP'))
        cache = cache_POP.out.cache.first()
        reference_archive = cache_POP.out.archive
    }
    results = reference_archive
    checkpoints = reference_archive
    if (params.steps.ascat) {
        ascat(ctx, stageScripts(scripts, 'ascat'), cache, [])
        checkpoints = checkpoints.concat(ascat.out.archive)
        results = results.concat(ascat.out.archive)
    }
    if (params.steps.pindel) {
        if (params.dependency_optimize) {
            input_roles = tumour_bas.map { tuple(1,it) }.mix(normal_bas.map { tuple(2,it) })
            PindelInputEarly(raw_ctx, stageScripts(scripts, 'PindelInputEarly'), cache, input_roles)
            input_parts = PindelInputEarly.out.parts
        } else {
            PindelInput(ctx, stageScripts(scripts, 'PindelInput'), cache, Channel.of(1,2))
            input_parts = PindelInput.out.parts
        }
        pindel_inputs = input_parts.toSortedList { a,b -> a[0] <=> b[0] }.map { rows -> rows.collect { it[1] } }
        PindelPlan(ctx, stageScripts(scripts, 'PindelPlan'), cache, pindel_inputs)
        pindel_plan = PindelPlan.out.plan.first()
        pindel_jobs = pindel_plan.flatMap { f ->
            def plan = new groovy.json.JsonSlurper().parseText(f.text)
            (1..(plan.jobs as int)).collect { i -> tuple(i, plan.chromosomes[i-1]) }
        }
        PindelChromosome(ctx, stageScripts(scripts, 'PindelChromosome'), cache, pindel_inputs, pindel_plan, pindel_jobs)
        if (params.dependency_optimize) {
            PindelFinalize(ctx, stageScripts(scripts, 'PindelFinalize'), cache, pindel_plan, PindelChromosome.out.parts.map { i,p -> p }.collect())
            PindelPack(PindelFinalize.out.result)
            pindel_archive = PindelPack.out.archive
            pindel_germline = PindelFinalize.out.germline
            checkpoints = checkpoints.concat(pindel_germline)
        } else {
            PindelGather(ctx, stageScripts(scripts, 'PindelGather'), cache, pindel_plan, PindelChromosome.out.parts.map { i,p -> p }.collect())
            pindel_archive = PindelGather.out.archive
            pindel_germline = PindelGather.out.archive
        }
        checkpoints = checkpoints.concat(pindel_archive)
        pindel_final = pindel_archive
        if (!params.skipannot) {
            cgpPindel_annot(ctx, stageScripts(scripts, 'cgpPindel_annot'), cache, pindel_archive)
            checkpoints = checkpoints.concat(cgpPindel_annot.out.archive)
            pindel_final = cgpPindel_annot.out.archive
        }
        results = results.concat(pindel_final)
    }
    if (params.steps.caveman) {
        CaVEMan_setup(ctx, stageScripts(scripts, 'CaVEMan_setup'), cache, [])
        checkpoints = checkpoints.concat(CaVEMan_setup.out.archive)
        CaVEMan_split(ctx, stageScripts(scripts, 'CaVEMan_split'), cache, CaVEMan_setup.out.archive)
        checkpoints = checkpoints.concat(CaVEMan_split.out.archive)
        CaVEMan_plan(ctx, stageScripts(scripts, 'CaVEMan_plan'), cache, CaVEMan_split.out.archive.concat(ascat.out.archive).collect())
        scatter_context = CaVEMan_plan.out.context.first()
        scatter_plan = CaVEMan_plan.out.plan.map { f -> new groovy.json.JsonSlurper().parseText(f.text) }.first()
        mstep_jobs = scatter_context.combine(scatter_plan).flatMap { context, plan ->
            (1..(plan.mstep_groups as int)).collect { index -> tuple(index, context) }
        }
        CaVEMan_Mstep(ctx, stageScripts(scripts, 'CaVEMan_Mstep'), cache, mstep_jobs)
        CaVEMan_ModelMerge(ctx, stageScripts(scripts, 'CaVEMan_ModelMerge'), cache, scatter_context, CaVEMan_Mstep.out.parts.map { i,p -> p }.collect())
        model = CaVEMan_ModelMerge.out.model.first()
        estep_jobs = scatter_context.combine(scatter_plan).flatMap { context, plan ->
            (1..(plan.estep_groups as int)).collect { index -> tuple(index, context) }
        }
        CaVEMan_Estep(ctx, stageScripts(scripts, 'CaVEMan_Estep'), cache, estep_jobs, model)
        CaVEMan_Gather(ctx, stageScripts(scripts, 'CaVEMan_Gather'), cache, scatter_context, CaVEMan_Estep.out.parts.map { i,p -> p }.collect(), model)
        checkpoints = checkpoints.concat(CaVEMan_plan.out.context,
            CaVEMan_Mstep.out.parts.map { i,p -> p }, CaVEMan_ModelMerge.out.model,
            CaVEMan_Estep.out.parts.map { i,p -> p }, CaVEMan_Gather.out.archive)
        flag_support = ascat.out.archive.concat(pindel_germline).collect()
        CaVEManFlagPlan(ctx, stageScripts(scripts, 'CaVEManFlagPlan'), cache, CaVEMan_Gather.out.archive.concat(ascat.out.archive, pindel_germline).collect())
        flag_jobs = CaVEManFlagPlan.out.splits.flatten().map { p -> tuple(p.name.replace('flag-group-', '').replace('.tar.gz', '').toInteger(), p) }
        CaVEManFlagGroup(ctx, stageScripts(scripts, 'CaVEManFlagGroup'), cache, flag_support, CaVEManFlagPlan.out.options.first(), flag_jobs)
        CaVEMan_flag(ctx, stageScripts(scripts, 'CaVEMan_flag'), cache, CaVEManFlagPlan.out.archive.concat(CaVEManFlagGroup.out.archive, ascat.out.archive, pindel_germline).collect())
        checkpoints = checkpoints.concat(CaVEMan_flag.out.archive)
        caveman_final = CaVEMan_flag.out.archive
        if (!params.skipannot) {
            CaVEMan_annot(ctx, stageScripts(scripts, 'CaVEMan_annot'), cache, CaVEMan_flag.out.archive)
            checkpoints = checkpoints.concat(CaVEMan_annot.out.archive)
            caveman_final = CaVEMan_annot.out.archive
        }
        results = results.concat(caveman_final)
    }
    if (params.steps.brass) {
        // Independent task directories; the driver rejects conflicting branch files.
        BRASSInput30(Channel.of(1,2), ctx, stageScripts(scripts, 'BRASSInput30'), cache, [])
        BRASSCoverage(ctx, stageScripts(scripts, 'BRASSCoverage'), cache, [])
        BRASSInsertSize(ctx, stageScripts(scripts, 'BRASSInsertSize'), cache, [])
        BRASSGroup(ctx, stageScripts(scripts, 'BRASSGroup'), cache, BRASSInput30.out.parts.map { i,p -> p }.collect())
        BRASSNormCN(ctx, stageScripts(scripts, 'BRASSNormCN'), cache, BRASSCoverage.out.archive.concat(ascat.out.archive).collect())
        BRASSPlan(ctx, stageScripts(scripts, 'BRASSPlan'), cache,
            BRASSGroup.out.archive.concat(BRASSInsertSize.out.archive, BRASSNormCN.out.archive, ascat.out.archive).collect())
        brass_seed = BRASSPlan.out.seed.first()
        brass_plan = BRASSPlan.out.plan.map { f -> new groovy.json.JsonSlurper().parseText(f.text) }.first()
        brass_jobs = brass_seed.combine(brass_plan).flatMap { seed, plan ->
            (1..plan.groups.size()).collect { i -> tuple(i, seed) }
        }
        BRASSAssemble(ctx, stageScripts(scripts, 'BRASSAssemble'), cache, brass_jobs)
        BRASSGather(ctx, stageScripts(scripts, 'BRASSGather'), cache,
            BRASSAssemble.out.parts.map { i,p -> p }.concat(BRASSPlan.out.archive, ascat.out.archive).collect())
        checkpoints = checkpoints.concat(BRASSInput30.out.parts.map { i,p -> p }, BRASSCoverage.out.archive,
            BRASSInsertSize.out.archive, BRASSGroup.out.archive, BRASSNormCN.out.archive,
            BRASSPlan.out.archive, BRASSPlan.out.seed, BRASSPlan.out.plan,
            BRASSAssemble.out.parts.map { i,p -> p }, BRASSGather.out.archive)
        results = results.concat(BRASSGather.out.archive)
    }
    if (params.steps.qc && !params.skipqc) {
        geno(ctx, stageScripts(scripts, 'geno'), cache, [])
        checkpoints = checkpoints.concat(geno.out.archive)
        verify_WT(ctx, stageScripts(scripts, 'verify_WT'), cache, [])
        checkpoints = checkpoints.concat(verify_WT.out.archive)
        verify_MT(ctx, stageScripts(scripts, 'verify_MT'), cache, ascat.out.archive)
        checkpoints = checkpoints.concat(verify_MT.out.archive)
        results = results.concat(geno.out.archive, verify_WT.out.archive, verify_MT.out.archive)
    }
    assemble(ctx, stageScripts(scripts, 'assemble'), results.collect())
    emit:
    intermediates = checkpoints.collect()
    run_params = assemble.out.run_params
    result_archive = assemble.out.result_archive
    timings = assemble.out.timings
    global_time = assemble.out.global_time
}
