nextflow.enable.dsl=2
include { stageScripts; BRASSInput30; BRASSCoverage; BRASSInsertSize; BRASSGroup; BRASSNormCN; BRASSPlan; BRASSAssemble; BRASSGather } from '../../modules/sanger-split'
workflow {
    ctx = Channel.value(tuple(*(1..9).collect { file("${projectDir}/fixtures/input${it}") }))
    scripts = file("${projectDir}/../../scripts/sanger-split")
    cache = Channel.value(file("${projectDir}/fixtures/cache"))
    ascat_archive = Channel.of(file("${projectDir}/fixtures/input1"))
        // Independent task directories; the driver rejects conflicting branch files.
        BRASSInput30(Channel.of(1,2), ctx, stageScripts(scripts, 'BRASSInput30'), cache, [])
        BRASSCoverage(ctx, stageScripts(scripts, 'BRASSCoverage'), cache, [])
        BRASSInsertSize(ctx, stageScripts(scripts, 'BRASSInsertSize'), cache, [])
        BRASSGroup(ctx, stageScripts(scripts, 'BRASSGroup'), cache, BRASSInput30.out.parts.map { i,p -> p }.collect())
        BRASSNormCN(ctx, stageScripts(scripts, 'BRASSNormCN'), cache, BRASSCoverage.out.archive.concat(ascat_archive).collect())
        BRASSPlan(ctx, stageScripts(scripts, 'BRASSPlan'), cache,
            BRASSGroup.out.archive.concat(BRASSInsertSize.out.archive, BRASSNormCN.out.archive, ascat_archive).collect())
        brass_seed = BRASSPlan.out.seed.first()
        brass_plan = BRASSPlan.out.plan.map { f -> new groovy.json.JsonSlurper().parseText(f.text) }.first()
        brass_jobs = brass_seed.combine(brass_plan).flatMap { seed, plan ->
            (1..plan.groups.size()).collect { i -> tuple(i, seed) }
        }
        BRASSAssemble(ctx, stageScripts(scripts, 'BRASSAssemble'), cache, brass_jobs)
        BRASSGather(ctx, stageScripts(scripts, 'BRASSGather'), cache,
            BRASSAssemble.out.parts.map { i,p -> p }.concat(BRASSPlan.out.archive, ascat_archive).collect())
        BRASSGather.out.archive.view()
    }
