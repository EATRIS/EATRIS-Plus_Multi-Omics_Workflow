# EATRIS-Plus Multi-omics Analysis Workflow
Analysis workflow used to analyze the cohort of healthy blood donors


# Prerequisites

### Multi-omics data set from Zenodo
- The data is stored in a Multi_Assay_Experiment object, which is used as input for the workflow
- The object can be downloaded here: https://doi.org/10.5281/zenodo.17514796


###  Install Nextflow using conda

 **Create a Conda Environment:** 
```bash
   conda create -n nextflow-env
   conda activate nextflow-env
   conda install -c bioconda nextflow
```

See also https://anaconda.org/bioconda/nextflow and https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.

### Singularity

For detailed Singularity installation instructions, please refer to the official [Singularity installation guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).


### Software containers
All software containers used in this workflow [can be obtained here](https://github.com/Xomics/Docker_containers) 


# Execute Analysis Workflow

For more instructions on how to run the workflow:
```
nextflow run main.nf --help
```

The typical command to run the workflow is:
```
nextflow run main.nf 
 	-c dre.config 
	--config_file config.yml
	--references references.bib 
	--mae_object /dir/MAE_object
	--container_dir /dir/containers/
	--output dir/of/choice
	
```
