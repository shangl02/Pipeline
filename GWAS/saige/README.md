# Nextflow pipeline to run SAIGE
This pipeline is for running GWAS association analysis using saige.
Saige has two steps, the first step is to fit a null model, the second step is to run the association analysis.

This folder has 4 pipelines.
1. **p01_main_chip.nf** single variant chip data GWAS.
2. **p02_main_wes.nf** single variant WES data GWAS.
3. **p03_main_wes_gene_saige45.nf** gene based WES data GWAS. (Use saige version 0.45)
4. **p04_main_wes_gene.nf** gene based WES data GWAS. (Use saige version > 1)

# Analysis for single variant GWAS
## Software container
We use singularity/docker image to run the analysis. Command to download the image for singularity is:

    singularity pull docker//shl198/saige_nextflow:v1.1.6

## Run pipeline
so we have two main file here for each file. The command to run the pipeline is:

    singularity run -B /:/media saige_nextflow.sif nextflow run /media/path/to/p01_main_wes.nf -c /media/path/to/p01_saige_wes.config -resume

# Analysis for gene based GWAS
For some reason, saige version > 1 reports error when running gene based analysis for UKBB data. So we use version 0.45 to run gene based GWAS.

## Software container
We use singularity/docker image to run the analysis. Command to download the image for singularity is:

    singularity pull docker//shl198/saige_nextflow:0.45

## Run pipeline
	singularity run -B /:/media saige_nextflow.sif nextflow run /media/path/to/p03_main_wes_gene_saige45.nf -c /media/path/to/p03_saige_wes_gene45.config -resume