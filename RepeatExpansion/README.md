# This folder has pipeline to detect repeat length using longread sequencing technology

This pipeline is used to build pipeline to detect repeat expansion for Nanopore Oxford long read sequencing. 

## Download docker container

	singularity pull docker://shl198/nanopore_re:v1

## How to run the pipeline
1 Define parameters in file p01_runNanoporePPL_Parameters.yaml <br/>
2 run the following command:
	
	singularity run -B /:/media \
		python /media/path/to/p01_runNanoporePPL.py \
				/media/path/to/p01_runNanoporePPL_Parameters.yaml