#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


process RUN_CV {

    publishDir "${params.output}/CV", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path references
    
    output:
    path 'cv_plots.html'

    """
    cp -L $project_dir/bin/features_overview.Rmd features_overview.Rmd 
    Rscript -e  "rmarkdown::render('features_overview.Rmd ', output_format = 'html_document', output_file = 'cv_plots.html', params = list(mae_hdf5_dir_path = '${mae_path}', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}


process LINEAR_MODELS {

    publishDir "${params.output}/RLM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path references
    
    output:
    path 'linear_models.html'

    """
    cp -L $project_dir/bin/rlm_processed_data.Rmd rlm_processed_data.Rmd
    Rscript -e  "rmarkdown::render('rlm_processed_data.Rmd', output_format = 'html_document', output_file = 'linear_models.html', params = list(mae_hdf5_dir_path = '${mae_path}', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}

process RUN_PCA {

    publishDir "${params.output}/PCA", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path references
    
    output:
    path 'pca_results.html'

    """
    cp -L $project_dir/bin/run_PCA_MAE.Rmd run_PCA_MAE.Rmd
    Rscript -e  "rmarkdown::render(run_PCA_MAE.Rmd', output_format = 'html_document', output_file = 'pca_results.html', params = list(mae_hdf5_dir_path = '${mae_path}', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}


process CELLTYPE_LEVELS {

    publishDir "${params.output}/Suppl", mode: 'copy', overwrite: true

    input:
    path mae_path
    path config_file
    path references
    
    output:
    path 'celltypes_results.html'

    """
    cp -L $project_dir/bin/celltype_levels.Rmd celltype_levels.Rmd
    Rscript -e  "rmarkdown::render(celltype_levels.Rmd', output_format = 'html_document', output_file = 'celltypes_results.html', params = list(mae_hdf5_dir_path = '${mae_path}', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}

process SINGLE_FEATURES {

    publishDir "${params.output}/Suppl", mode: 'copy', overwrite: true

    input:
    path mae_path
    path config_file
    path references
    
    output:
    path 'single_features_results.html'

    """
    cp -L $project_dir/bin/analyze_single_features.Rmd analyze_single_features.Rmd
    Rscript -e  "rmarkdown::render(analyze_single_features.Rmd', output_format = 'html_document', output_file = 'single_features_results.html', params = list(mae_hdf5_dir_path = '${mae_path}', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}




