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
	    	--mudata dir/to/mudata.hd5		
			--r_config path/to/r_config.yml
			--references path/to/references.bib
       Optional arguments:
	--container_dir				The directory where the required Singularity (.sif files) images are stored
	--sil_cutoff				Cutoff silhouette score used for filtering fused networks (default = 0.05)
	--min_n_cutoff				Minimal number of data modalities combinations to run SNF on (default = 2)
	--max_n_cutoff				Maximal number of data modalities combinations to run SNF on (default = 10)
        """
}




////////////////////////////////////////////////////
/* --            Input data files              -- */
////////////////////////////////////////////////////+


mae_object = Channel.fromPath("${params.mae_object}")
mudata = Channel.fromPath("${params.mudata}")  
temp_file = Channel.fromPath("${params.temp_file}")
r_config = Channel.fromPath("${params.config_file}")
references = Channel.fromPath("${params.references}")
rna_snp_file = Channel.fromPath("${params.rna_snp}")
snp_cpg_files = Channel.fromPath("${params.snp_cog_dir}/*")
r_config = Channel.fromPath("${params.config_file}")

////////////////////////////////////////////////////
/* --                  Modules                 -- */
////////////////////////////////////////////////////+


include { CELLTYPE_LEVELS; SINGLE_FEATURES } from './modules/single_omics.nf'
include { LINEAR_MODELS; RUN_CV; RUN_PCA } from './modules/single_omics.nf'
include { COMPARE_LIPIDOMICS; COMPARE_PROTEOMICS; SAMPLES_PER_OMICS } from './modules/single_omics.nf'
include { RUN_MOFA; RIBOSOMAL_RNA_RANKS } from './modules/mofa.nf'
include { RUN_PLS2; RUN_SNP_CPG } from './modules/pairwise_omics.nf'
include { PREPARE_MUDATA } from './modules/mudata.nf'
include { SINGLE_OMICS_CLUSTERING ; OMICS_COMBINATIONS ; SNF ; CONCAT_SIL_SCORES ; PLOT_SIL_SCORES; SNF_ANALYSIS; SNF_VENN_DIAGRAM } from './modules/snf'



////////////////////////////////////////////////////
/* --                 Workflow                 -- */
////////////////////////////////////////////////////+


workflow {

	// Show help message
	if (params.help) {
    	helpMessage()
    	exit 0
	}


	/// Read mudata object and filter samples/features
	//PREPARE_MUDATA(mudata)
	

	// General analysis 
	//CELLTYPE_LEVELS(params.mae_object, params.r_config, params.references)
	//SAMPLES_PER_OMICS(params.mae_object, params.r_config) 


	// Single omics analyses
	//RUN_CV(params.mae_object, params.r_config, params.references)
	LINEAR_MODELS(params.mae_object, params.r_config, params.references)
	//SINGLE_FEATURES(params.mae_object, params.r_config, params.references, LINEAR_MODELS.out[1])
	//RUN_PCA(params.mae_object, params.r_config, params.references)
	COMPARE_PROTEOMICS(params.mae_object, params.r_config, LINEAR_MODELS.out[1])
	COMPARE_LIPIDOMICS(params.mae_object, params.r_config, LINEAR_MODELS.out[1])


	// Multi omics analyses
	//RUN_MOFA(params.mae_object, params.r_config, params.references)
	//RIBOSOMAL_RNA_RANKS(params.mae_object, params.r_config, params.references)
	//RUN_PLS2(params.mae_object, params.r_config, params.references)
	RUN_SNP_CPG(params.mae_object, params.rna_snp_file, params.snp_cpg_files)

	// Similarity Network Fusion
	//SINGLE_OMICS_CLUSTERING(mudata)
	//OMICS_COMBINATIONS(mudata, params.min_n_combi, params.max_n_combi)
	//SNF(mudata.toList(), OMICS_COMBINATIONS.out.flatten(), params.sil_cutoff)
	//CONCAT_SIL_SCORES(SNF.out[0].toList())
	//PLOT_SIL_SCORES(CONCAT_SIL_SCORES.out)
	//SNF_ANALYSIS(SNF.out[1].flatten(), mudata.toList())
}
