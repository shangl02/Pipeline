# Nextflow pipeline to run SAIGE(version > 1)
This pipeline is for running GWAS association analysis using saige.
Saige has two steps, the first step is to fit a null model, the second step is to run the association analysis.

## Software container
We use singularity/docker image to run the analysis. Command to download the image for singularity is:

    singularity pull docker//shl198/saige_nextflow:v1.1.6

## Run pipeline
For Chip data, all of the variants for step1 are in one plink file. For WES data, variants for step1 are split by chromosome, so we have two main file here for each file. The command to run the pipeline is:

    singularity run -B /:/media saige_nextflow.sif nextflow run /media/path/to/main_wes.nf -c /media/path/to/saige_wes.config

