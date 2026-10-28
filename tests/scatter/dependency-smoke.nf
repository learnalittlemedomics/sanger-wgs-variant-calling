nextflow.enable.dsl=2
include { stageScripts; prepareEarly; cacheEarly; PindelInputEarly; PindelPlan; PindelChromosome; PindelFinalize; PindelPack } from '../../modules/sanger-split'
process BasMock {
 input:
 val role
 output:
 tuple val(role), path("input${role == 1 ? 4 : 7}.bas")
 script:
 """
 sleep ${role == 1 ? 1 : 8}
 echo bas > input${role == 1 ? 4 : 7}.bas
 """
}
workflow {
 f = (1..9).collect { file("${projectDir}/fixtures/input${it}") }
 scripts = file("${projectDir}/../../scripts/sanger-split")
 BasMock(Channel.of(1,2))
 roles = BasMock.out
 prepareEarly(f[0],f[1],f[2],f[5],f[8],f[3],f[4],f[6],f[7],stageScripts(scripts,'prepareEarly'))
 raw = prepareEarly.out.context.first()
 cacheEarly(raw,stageScripts(scripts,'cacheEarly'))
 cache=cacheEarly.out.cache.first()
 PindelInputEarly(raw,stageScripts(scripts,'PindelInputEarly'),cache,roles)
 inputs=PindelInputEarly.out.parts.toSortedList { a,b -> a[0] <=> b[0] }.map { rows -> rows.collect { it[1] } }
 tb=roles.filter { i,p -> i == 1 }.map { i,p -> p }
 nb=roles.filter { i,p -> i == 2 }.map { i,p -> p }
 ctx=raw.combine(tb).combine(nb).map { r,p,f,t,ti,n,ni,tb,nb -> tuple(r,p,f,t,ti,tb,n,ni,nb) }.first()
 PindelPlan(ctx,stageScripts(scripts,'PindelPlan'),cache,inputs)
 plan=PindelPlan.out.plan.first()
 jobs=plan.flatMap { p -> def q=new groovy.json.JsonSlurper().parseText(p.text); (1..q.jobs).collect { i -> tuple(i,q.chromosomes[i-1]) } }
 PindelChromosome(ctx,stageScripts(scripts,'PindelChromosome'),cache,inputs,plan,jobs)
 PindelFinalize(ctx,stageScripts(scripts,'PindelFinalize'),cache,plan,PindelChromosome.out.parts.map { i,p -> p }.collect())
 PindelPack(PindelFinalize.out.result)
 PindelFinalize.out.germline.view { "BED_READY: ${it}" }
}
