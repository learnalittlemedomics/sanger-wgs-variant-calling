nextflow.enable.dsl=2
params.publish_dir="${launchDir}/sanger-subset-output"
include { payload as pGenVarSnv; payload as pGenVarIndel; payload as pGenVarCnv; payload as pGenVarSv; payload as pGenVarSupp; payload as pGenQc } from './fixture'
include { verifiedPublish } from '../../modules/verified-publish'
workflow {
local_mode=true
sanger_steps=[ascat:true,pindel:false,caveman:false,brass:false]
pGenVarSnv('SangerWgs_pGenVarSnv')
pGenVarIndel('SangerWgs_pGenVarIndel')
pGenVarCnv('SangerWgs_pGenVarCnv')
pGenVarSv('SangerWgs_pGenVarSv')
pGenVarSupp('SangerWgs_pGenVarSupp')
pGenQc('SangerWgs_pGenQc')
        // Completion and cleanup depend on durable, verified publication.
        if (local_mode) {
            publication_groups = pGenVarSupp.out.payload.map { [role: 'SangerWgs_pGenVarSupp', payload: it.toString()] }
            publication_groups = publication_groups.concat(pGenQc.out.payload.map { [role: 'SangerWgs_pGenQc', payload: it.toString()] })
            publication_expected = ['SangerWgs_pGenVarSupp', 'SangerWgs_pGenQc']
            if (sanger_steps.ascat) {
                publication_groups = publication_groups.concat(pGenVarCnv.out.payload.map { [role: 'SangerWgs_pGenVarCnv', payload: it.toString()] })
                publication_expected = publication_expected + ['SangerWgs_pGenVarCnv']
            }
            if (sanger_steps.pindel) {
                publication_groups = publication_groups.concat(pGenVarIndel.out.payload.map { [role: 'SangerWgs_pGenVarIndel', payload: it.toString()] })
                publication_expected = publication_expected + ['SangerWgs_pGenVarIndel']
            }
            if (sanger_steps.caveman) {
                publication_groups = publication_groups.concat(pGenVarSnv.out.payload.map { [role: 'SangerWgs_pGenVarSnv', payload: it.toString()] })
                publication_expected = publication_expected + ['SangerWgs_pGenVarSnv']
            }
            if (sanger_steps.brass) {
                publication_groups = publication_groups.concat(pGenVarSv.out.payload.map { [role: 'SangerWgs_pGenVarSv', payload: it.toString()] })
                publication_expected = publication_expected + ['SangerWgs_pGenVarSv']
            }
            verifiedPublish(publication_groups.collect(flat: false), publication_expected, params.publish_dir)
        }

}
