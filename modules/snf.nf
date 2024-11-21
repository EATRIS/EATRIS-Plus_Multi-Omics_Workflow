#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


process SINGLE_OMICS_CLUSTERING {

    publishDir "${params.output}/SNF", mode: 'copy', overwrite: true

    input:
    path mudata_object
    
    output:
    path 'single_omics_clustering_out.ipynb'

    """
    python -m ipykernel install --name base --user
    papermill $project_dir/bin/single_omics_clustering.ipynb single_omics_clustering_out.ipynb -p mudata_object ${mudata_object}  -k base
    """
    
}


process OMICS_COMBINATIONS {

    input:
    path mudata_object
    val min_n_combi
    val max_n_combi
    
    output:
    path '*.txt'

    """
    python $project_dir/bin/find_omics_combinations.py ${mudata_object} ${min_n_combi} ${max_n_combi}
    """
    
}


process SNF {

	label 'full_resources'

	input:
	path mudata_object
	path omics_names
	val sil_cutoff

	output:
	path '*_sil_scores.tsv'
	path '*_fused_network.csv', emit: test_channel optional true
    
	"""
	python $project_dir/bin/snf_run.py  ${mudata_object} ${omics_names} ${sil_cutoff} 
	"""
}

process CONCAT_SIL_SCORES {

    publishDir "${params.output}/SNF", mode: 'copy', overwrite: true

	label 'full_resources'

	input:
	path chunks
    
	output:
	path 'sil_scores_combined.csv'
    
	"""
	cat ${chunks} > sil_scores_combined.csv
	"""
}

process PLOT_SIL_SCORES {

	publishDir "${params.output}/SNF", mode: 'copy', overwrite: true

	label 'full_resources'

	input:	
	path sil_combined
    
	output:
	path 'get_top_fused_networks_out.ipynb'
    	
	""" 
	papermill $project_dir/bin/plot_silhouette_scores.ipynb  get_top_fused_networks_out.ipynb  -p sil_combined_path ${sil_combined}   -k base
	""" 

}

process SNF_ANALYSIS {

	label 'full_resources'

	publishDir "${params.output}/SNF/Analysis", mode: 'copy', overwrite: true

	input:
	path array
	path mudata_object

	output:
	path '*.ipynb'
    path 'labels/SNF_labels_*.csv'

	"""
	python -m ipykernel install --name base --user
	papermill $project_dir/bin/snf_analysis.ipynb  ${array}.ipynb  -p mudata_object ${mudata_object} -p snf_df_path ${array} --stdout-file ${array}_test.ipynb  -k base
	"""

}

process SNF_VENN_DIAGRAM {

    publishDir "${params.output}/SNF/VennDiagrams", mode: 'copy', overwrite: true

	input:
	path labels_dataframes

	output:
	path '*.ipynb'


	"""
	python -m ipykernel install --name base --user
	papermill $project_dir/bin/snf_venn_diagram.ipynb  venn_diagram.ipynb -p labels_paths ${labels_dataframes} -k base
	"""


}
