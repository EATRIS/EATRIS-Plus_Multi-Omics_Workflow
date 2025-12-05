#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


process RUN_CV {

    publishDir "${params.output}/CV", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path references
    
    output:
    path 'cv_plots.html'

    """
    cp -L $project_dir/bin/features_overview.Rmd features_overview.Rmd 
    Rscript -e  "rmarkdown::render('features_overview.Rmd', output_format = 'html_document', output_file = 'cv_plots.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}


process LINEAR_MODELS {

    publishDir "${params.output}/RLM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path references
    
    output:
    path 'linear_models.html'
    path 'rlm_signif_hits.csv'

    """
    cp -L $project_dir/bin/rlm_processed_data.Rmd rlm_processed_data.Rmd
    Rscript -e  "rmarkdown::render('rlm_processed_data.Rmd', output_format = 'html_document', output_file = 'linear_models.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}

process RUN_PCA {

    publishDir "${params.output}/PCA", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path references
    
    output:
    path 'pca_results.html'

    """
    cp -L $project_dir/bin/run_PCA_MAE.Rmd run_PCA_MAE.Rmd
    Rscript -e  "rmarkdown::render('run_PCA_MAE.Rmd', output_format = 'html_document', output_file = 'pca_results.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}


process CELLTYPE_LEVELS {

    publishDir "${params.output}/Suppl", mode: 'copy', overwrite: true

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path references
    
    output:
    path 'celltypes_results.html'

    """
    cp -L $project_dir/bin/celltype_levels.Rmd celltype_levels.Rmd
    cp -L $project_dir/SupportFunctions/ViolinPlot.R ViolinPlot.R
    Rscript -e  "rmarkdown::render('celltype_levels.Rmd', output_format = 'html_document', output_file = 'celltypes_results.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', references = '${references}'  ))"
    """
    
}

process SINGLE_FEATURES {

    publishDir "${params.output}/Suppl", mode: 'copy', overwrite: true

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path references
    path rlm_signif_hits
    
    output:
    path 'single_features_results.html'

    """
    cp -L $project_dir/bin/analyze_single_features.Rmd analyze_single_features.Rmd
    cp -L $project_dir/SupportFunctions/ViolinPlot.R ViolinPlot.R
    cp -L $project_dir/SupportFunctions/PlotCor2.R PlotCor2.R
    cp -L $project_dir/SupportFunctions/themes_and_colors.R themes_and_colors.R
    Rscript -e  "rmarkdown::render('analyze_single_features.Rmd', output_format = 'html_document', output_file = 'single_features_results.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', references = '${references}', rlm_sign_hits = '${rlm_signif_hits}'  ))"
    """
    
}


process COMPARE_PROTEOMICS {

    publishDir "${params.output}/PLATFORMS", mode: 'copy', overwrite: true

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path rlm_signif_hits
    
    output:
    path 'compare_prot_platforms.html'

    """
    cp -L $project_dir/bin/compare_proteomics_platforms.Rmd compare_proteomics_platforms.Rmd
    Rscript -e  "rmarkdown::render('compare_proteomics_platforms.Rmd', output_format = 'html_document', output_file = 'compare_prot_platforms.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', rlm_sign_hits = '${rlm_signif_hits}'  ))"
    """
}

process COMPARE_LIPIDOMICS {

    publishDir "${params.output}/PLATFORMS", mode: 'copy', overwrite: true

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    path rlm_signif_hits
    
    output:
    path 'compare_lip_platforms.html'

    """
    cp -L $project_dir/bin/compare_lipidomics_platforms.Rmd compare_lipidomics_platforms.Rmd
    Rscript -e  "rmarkdown::render('compare_lipidomics_platforms.Rmd', output_format = 'html_document', output_file = 'compare_lip_platforms.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}', rlm_sign_hits = '${rlm_signif_hits}'  ))"
    """
}


process SAMPLES_PER_OMICS {

    publishDir "${params.output}/STUDY_DESIGN", mode: 'copy', overwrite: true

    input:
    path mae_path, stageAs: 'mae_hdf5'
    path config_file
    
    output:
    path 'study_design.html'
    path 'plot_samples_per_omics.png'

    """
    cp -L $project_dir/bin/plot_samples_per_omics.Rmd plot_samples_per_omics.Rmd
    Rscript -e  "rmarkdown::render('plot_samples_per_omics.Rmd', output_format = 'html_document', output_file = 'study_design.html', params = list(mae_hdf5_dir_path = 'mae_hdf5', config_file = '${config_file}'))"
    """
    
}




