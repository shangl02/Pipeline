GATK pipeline for whole genome/exome sequencing
========================================
It's wirtten in nextflow.
There are two parts for this pipeline.
1. Generate gvcf files for all the samples.
	- a) define all the parameters in file p01_GATK_get_gvcf_Parameters.config.
	- b) in terminal run **nextflow run p01_GATK_get_gvcf.nf -c p01_GATK_get_gvcf_Parameters.config -resume**

2. Do VQSR for the mutations.
	- a) define all the parameters in file p02_GATK_VQSR_Parameters.config.
	- b)**nextflow run p02_GATK_VQSR.nf -c p02_GATK_VQSR_Parameters.config -resume**

3. Samples with multiple lanes.
	- Files starts with p03 perform the same function as files that starts with p01. It's used for samples with multiple lanes. As suggested by GATK https://gatkforums.broadinstitute.org/gatk/discussion/3060/how-should-i-pre-process-data-from-multiplexed-sequencing-and-multi-library-designs. First align the raw fastq files from each lane separately with differetn read groups, then merge the lanes that belong to the same group together.

4. Identify DeNovo variants using GATK
	- This pipeline follows the suggested pipeline here https://gatk.broadinstitute.org/hc/en-us/articles/360035531432?id=11074 to find the DeNovo mutations.
	- How to run it:
		- **nextflow run p04_GATK_DeNovo.nf -c p04_GATK_DeNovo.config -resume**

Softwares that were installed in the docker file
------------------------------------
* GATK 4.1.4.1
* fastqc 0.11.9
* bwa 0.7.17
* samtools 1.9
* tabix 0.2.5
* trimmomatic 0.39
* bcftools 1.8
* picardtools 2.21.2
