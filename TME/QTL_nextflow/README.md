# This pipeline aims to run the QTL analaysis for identifying the association between variants and immune cell abundance. It assumes the genotype already been imputed for genotype array/WES, or the genotype is from WGS/WES and doesn't need to impute.

## Input files:

Before running the pipeline, you need to prepare some files.
* Vcf file in gzip format that has genotypes of interested cohort.
* A tabix index gzipped file that have annotation information of the variants in vcf file, annotation can be done using VEP, please check [here](https://github.com/lis262/Variant_Annotation) to generate the annotation file.
* An id mapping file that has the follwing columns. "vcfid, sid, cytoid,Age,Sex". vcfid is the sample name in vcf file, sid is subject id, cytoid is the sample name in cell abundance file. 
* Cell abunance file with rows as samples, columns as cell types.
* Reference plink bfile that will be used to predict ethnicity and calculate Pinciple component for covariates in the model. You can use 1000 genomics file, or use the bfile coverted from the vcf file.
* high LD region file which can be downloaded [here](https://dougspeed.com/high-ld-regions/), it's for hg19, if you want to use hg38, you need to liftover it. I provide both hg19 and hg38 versions here in the repository.


## Docker contianer download

	singularity pull docker://shl198/tme:v1


## How to run the pipeline
1 Define all parameters in tme_ppl.config
2 Run the following comand:

	singularity run -B /:/media tem_v1.sif \
		nextflow /media/path/to/tme_ppl.nf \
			-c /media/path/to/tme_ppl.config