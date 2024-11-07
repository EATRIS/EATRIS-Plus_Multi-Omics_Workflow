# EATRIS-Plus Multi-omics Analysis Workflow
Analysis workflow used to analyze the cohort of healthy blood donors


# Requirements

Nextflow and Singularity are required to be installed. To install them on CentOS Linux 7 enable [EPEL](https://docs.fedoraproject.org/en-US/epel/) and run

```
sudo yum -y install conda.noarch
sudo yum -y install singularity.x86_64
```

## Installing Nextflow
```
conda create --name nf
conda activate nf
conda install -c bioconda nextflow:22.04.0
```

See also https://anaconda.org/bioconda/nextflow and https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.

## Run Workflow

For more instructions on how to run the workflow:
```
nextflow run main.nf --help
```

The typical command to run the workflow is:
```
nextflow run main.nf 
	--output dir/of/choice
 	-c dre.config 
	--config_file config.yml 
	--mae_object /dir/MAE_object
	--container_dir /dir/Singularity/containers

```


## 