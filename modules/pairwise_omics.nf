#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


process RUN_PLS2 {

    publishDir "${params.output}/PLS2", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path references
    
    output:
    path 'pls2.html'

    """
    cp -L $project_dir/bin/run_mixOmics_PLS2.Rmd run_mixOmics_PLS2.Rmd
    Rscript -e  "rmarkdown::render('run_mixOmics_PLS2.Rmd', output_format = 'html_document', output_file = 'pls2.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}



process RUN_SNP_CPG {

    publishDir "${params.output}/SNP_CPG", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path rna_file
    path snp_cpg_files
    
    output:
    path 'snp_cpg_analysis.html'

    """
    cp -L $project_dir/bin/analyze_SNP-CpGs.Rmd analyze_SNP-CpGs.Rmd
    Rscript -e  "rmarkdown::render('analyze_SNP-CpGs.Rmd', output_format = 'html_document', output_file = 'snp_cpg_analysis.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', rna_lm_results = '${rna_file}', SNP_CpG_Dir = '${snp_cpg_files}'  ))"
    """
    
}


