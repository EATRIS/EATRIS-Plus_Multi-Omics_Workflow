#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


process RUN_PLS2 {

    publishDir "${params.output}/PLS2", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path references
    
    output:
    path 'pls2.html'

    """
    cp -L $project_dir/bin/run_PLS2.Rmd run_PLS2.Rmd
    Rscript -e  "rmarkdown::render('run_PLS2.Rmd', output_format = 'html_document', output_file = 'pls2.html', params = list(mae_hdf5_dir_path = '${mae_path}', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}


