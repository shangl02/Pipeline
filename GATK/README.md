GATK pipeline for whole exome sequencing
========================================
It's wirtten in nextflow.
There are two parts for this pipeline.
1. Generate gvcf files for all the samples.
	a) define all the parameters in file p01_GATK_get_gvcf_Parameters.config.
	b) in terminal run **nextflow run p01_GATK_get_gvcf.nf -c p01_GATK_get_gvcf_Parameters.config -resume**

2. Do VQSR for the mutations.
	a) define all the parameters in file p02_GATK_VQSR_Parameters.config.
	b)**nextflow run p02_GATK_VQSR.nf -c p02_GATK_VQSR_Parameters.config -resume**


Softwares that needs to be installed
------------------------------------
* GATK
* bwa
* samtools
* tabix
* trimmomatic
* bcftools