#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process PREPARE_MUDATA {


	input: 
	path mudata_object

	output:
	path 'mudata_preprocessed.h5mu'

	"""
	python3 $project_dir/bin/prepare_mudata.py ${mudata_object} mudata_preprocessed.h5mu 
	"""
}
