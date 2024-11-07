#!/usr/bin/env nextflow

nextflow.enable.dsl=2

project_dir = projectDir


////////////////////////////////////////////////////
/*    --               Functions               -- */
////////////////////////////////////////////////////+


def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf 
            --output dir/of/choice
            --mae_object path/to/mae_object

       Optional arguments:
       --container_dir				The directory where the required Singularity (.sif files) images are stored
        """
}




////////////////////////////////////////////////////
/* --            Input data files              -- */
////////////////////////////////////////////////////+


mae_object = Channel.fromPath("${params.mae_object}")
temp_file = Channel.fromPath("${params.temp_file}")
r_config = Channel.fromPath("${params.config_file}")
references = Channel.fromPath("${params.references}")

////////////////////////////////////////////////////
/* --                  Modules                 -- */
////////////////////////////////////////////////////+


include { RUN_CV } from './modules/single_omics.nf'
//include { LINEAR_MODELS } from './modules/single_omics.nf'
//include { RUN_PCA} from './modules/single_omics.nf'
//include { RUN_MOFA } from './modules/mofa.nf'
//include { RUN_PLS2 } from './modules/pls2.nf'



////////////////////////////////////////////////////
/* --                 Workflow                 -- */
////////////////////////////////////////////////////+


workflow {

	// Show help message
	if (params.help) {
    	helpMessage()
    	exit 0
	}


	// Single omics analyses
	RUN_CV(mae_object, r_config, references)
	//LINEAR_MODELS(params.mae_object, params.r_config, params.references)
	//RUN_PCA(params.mae_object, params.r_config, params.references)

	// Multi omics analyses
	//RUN_MOFA(params.mae_object, params.r_config, params.references)
	//RUN_PLS2(params.mae_object, params.r_config, params.references)
}
