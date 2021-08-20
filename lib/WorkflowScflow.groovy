//
// This file holds several functions specific to the workflow/scflow.nf in the nf-core/scflow pipeline
//
class WorkflowScflow {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        if (params.dge_mast_method == 'glmer' && params.dge_random_effects_var == 'NULL') {
            log.error "The glmer method requires a random effects variable."
            System.exit(1)
        }

        if (params.dge_mast_method == 'glm' && params.dge_random_effects_var != 'NULL') {
            log.error "The glm method can not fit a random effects variable."
            System.exit(1)
        }
    }
}
