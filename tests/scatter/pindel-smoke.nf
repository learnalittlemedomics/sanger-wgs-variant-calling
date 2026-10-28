nextflow.enable.dsl=2
include { stageScripts; PindelInput; PindelPlan; PindelChromosome; PindelGather } from '../../modules/sanger-split'
workflow {
 ctx = Channel.value(tuple(*(1..9).collect { file("${projectDir}/fixtures/input${it}") }))
 scripts = file("${projectDir}/../../scripts/sanger-split")
 cache = Channel.value(file("${projectDir}/fixtures/cache"))
 PindelInput(ctx, stageScripts(scripts,'PindelInput'), cache, Channel.of(1,2))
 inputs = PindelInput.out.parts.toSortedList { a,b -> a[0] <=> b[0] }.map { rows -> rows.collect { it[1] } }
 PindelPlan(ctx, stageScripts(scripts,'PindelPlan'), cache, inputs)
 plan = PindelPlan.out.plan.first()
 jobs = plan.flatMap { f -> def p = new groovy.json.JsonSlurper().parseText(f.text); (1..(p.jobs as int)).collect { i -> tuple(i,p.chromosomes[i-1]) } }
 PindelChromosome(ctx, stageScripts(scripts,'PindelChromosome'), cache, inputs, plan, jobs)
 PindelGather(ctx, stageScripts(scripts,'PindelGather'), cache, plan, PindelChromosome.out.parts.map { i,p -> p }.collect())
}
