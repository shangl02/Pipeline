# Pipeline to identify structure variants


## Download container using the following command

	singularity pull docker://shl198/smoove_nextflow:v0.2.5

The pipeline uses tool [smoove](https://github.com/brentp/smoove) which uses lumpy to call the structure variants.



## How  to run the pipeline
1 Define parameters in lumpy_smoove.config
2 Secondin terminal run command

	singularity run -B /:/media smoove_nextflow_v0.2.5.sif \
	  	nextflow run /media/path/to/lumpy_smoove.nf \
	  		-c /media/path/to/lumpy_smoove.config
