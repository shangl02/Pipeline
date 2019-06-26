GATK pipeline for whole exome sequencing
========================================
It's wirtten in nextflow.
The way to run the pipeline:
* 1. define all the parameters in GAK_Parameters.config.
* 2. in terminal type the following command: 
     nextflow run GATK_human.nf -c GATK_Parameters.config -resume

