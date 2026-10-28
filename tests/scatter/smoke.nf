nextflow.enable.dsl=2
// Exercise the actual production processes' staging, dynamic resources and retry.
// No CRAM, container, seff or Slurm is used with -stub-run and smoke.config.
include { stageScripts; CaVEMan_plan; CaVEMan_Mstep; CaVEMan_ModelMerge; CaVEMan_Estep; CaVEMan_Gather } from '../../modules/sanger-split' params(scatter: [target_hours:1.5, workers_per_job:6, worker_steps:[6,6,4,2], group_cpu_steps:[24,30,30,30], small_cpu_steps:[1,2,3,5], memory_gb_per_cpu:3.6, max_retries:3])
workflow {
    ctx = Channel.value(tuple(*(1..9).collect { file("${projectDir}/fixtures/input${it}") }))
    scripts = file("${projectDir}/../../scripts/sanger-split")
    cache = Channel.value(file("${projectDir}/fixtures/cache"))
    CaVEMan_plan(ctx, stageScripts(scripts, 'CaVEMan_plan'), cache, [])
    context = CaVEMan_plan.out.context.first()
    plan = CaVEMan_plan.out.plan.map { f -> new groovy.json.JsonSlurper().parseText(f.text) }.first()
    mjobs = context.combine(plan).flatMap { c,p -> (1..(p.mstep_groups as int)).collect { i -> tuple(i,c) } }
    CaVEMan_Mstep(ctx, stageScripts(scripts, 'CaVEMan_Mstep'), cache, mjobs)
    CaVEMan_ModelMerge(ctx, stageScripts(scripts, 'CaVEMan_ModelMerge'), cache, context, CaVEMan_Mstep.out.parts.map { i,p -> p }.collect())
    model = CaVEMan_ModelMerge.out.model.first()
    ejobs = context.combine(plan).flatMap { c,p -> (1..(p.estep_groups as int)).collect { i -> tuple(i,c) } }
    CaVEMan_Estep(ctx, stageScripts(scripts, 'CaVEMan_Estep'), cache, ejobs, model)
    CaVEMan_Gather(ctx, stageScripts(scripts, 'CaVEMan_Gather'), cache, context, CaVEMan_Estep.out.parts.map { i,p -> p }.collect(), model)
}
