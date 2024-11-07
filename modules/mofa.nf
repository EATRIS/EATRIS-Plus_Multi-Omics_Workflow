#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


process RUN_MOFA {

    publishDir "${params.output}/MOFA", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path references
    
    output:
    path 'mofa.html'

    """
    cp -L $project_dir/bin/run_MOFA.Rmd run_MOFA.Rmd
    Rscript -e  "rmarkdown::render('run_MOFA.Rmd', output_format = 'html_document', output_file = 'mofa.html', params = list(mae_hdf5_dir_path = '${mae_path}', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}


