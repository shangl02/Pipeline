This pipeline aims to run the QTL analaysis for identifying the association between variants and immune cell abundance. It assumes the genotype already been imputed for genotype array/WES, or the genotype is from WGS/WES and doesn't need to impute.

Input files:
------------

Before running the pipeline, you need to prepare some files.
* Vcf file in gzip format that hhas genotypes of interested cohort.
* A tabix index gzip file that have annotation information of the variants in vcf file.
* An id mapping file that has the follwing columns. "vcfid, sid, cytoid,Age,Sex". vcfid is the sample name in vcf file, sid is subject id, cytoid is the sample name in cell abundance file. 
* Cell abuance file with rows as samples, columns as cell types.
* Reference plink bfile that will be used to predict ethnicity and calculate Pinciple component for covariates in the model. You can use 1000 genomics file, or use the bfile coverted from the vcf file.


Softwares needed:
-----------------
* python3: some required packages are in file requirements.
* plink1: binary version.
* plink2: binary version.
* flashPCA: binary version.
* matrixQTL: will provide docker version in the future.


### Step 1) Transfer vcf to plink format and Filter the variants**

How to run the command:

    python m01_process_gt.py \
              -i vcf.gz \       # vcf file
              -p work_path \    # path to work on
              --af 0.01 \       # if don't specify, default is 0.01
              --plink2 path/to/plink2/binary \
              --plink1 path/to/plink1/binary

After this command, a folder named **plink** will be created in the work_path and plink files with prefix **germ_rmSNP38** will be created in it.

### Step 2) Clean samples

Thie step is to remove samples that come from the same person, we keep the one with lower SNP missing rate.

	python m02_clean_samples.py \
	         -p work_path \   # should be the same path as step1
	         --id_map id_mapping_file \
	         --plink2 path/to/plink2/binary \
             --plink1 path/to/plink1/binary

The id_mapping_file needs to have 6 columns, **vcfid, sid, cytoid, Sex, Age, Race.** sid is subject id. cytoid is the cytoreason sample id. This command creates a subfolder **f01_qc** in the plink folder created in step1. It remove the samples with duplicate WGS data and then it updates the sample name by changing the vcfid to DNA_sid, with plink bfile with prefix **germ_newSp38** as output. Finally it transfers the plink file to a table with rows as SNPS and columns as samples and values are genotype. The file name is **germ_newSp38.traw.gz**

### Step 3) Calculate PC forthe cohort
This step is used to calculate the principle component of SNPs to predict the ethnicity. It uses 1k genomeics as refrence.

	python m03_calculate_PC.py \
		-r ref_plink_bfile_pre \  # usually 1k genomics plink bfile prefix
		-p work_path \  # same as previous steps
		--plink2 path/to/plink2/binary \
        --plink1 path/to/plink1/binary \
        -f flashPCA     # path to flashPCA bianry tool

This step creats a folder nameed **f02_pca** in plink folder. It gets overlapped SNPs between cohort and 1k genomics and use these SNPs from 1k genomics to calculate PCA and then project to the cohorts to get the PC for the cohorts.

### Step 4) Prepare covarate file and cytoreason file for MatrixQTL anlaysis

This step will add 'DNA_' as prefix for subject id in id mapping file. Currently it conly considers first 10 PC, Age and Sex, if you want to include more covariates, you can modify the code.

	python m04_get_covariate_file.py \
		-p work_path \         # same as previous steps
		--id_map id_map_fn \   # file id map file used in step2
		-d deconvolute_cell_fn # file with cell type abundance data 
After the command, it will generate two files: one named 'covariates.txt' with rows as covariates and columns as samples, another named 'cyto.txt' with rows as cell types, columns as samples. Please check if the sample names are the same as in the file **germ_newSp38.traw.gz** in step 2.

### Step 5) Prepare files for MatrixQTL
Sometimes cytoreason file, covaraite file and sample file don't have the same samles. So this step will first get the overlapped sample names and then prepare the input files for matrixQTL.

	python m05_getinput4matrixQTL.py \	
		-c covaraite_fn \      # generated in step 4
		-d deconvolution_fn \  # deconvolution file
		-p work_path \         # same as previous steps

This step will generate 3 files: SNP.txt, cyto.txt, covariates.txt which will be the input for matrixQTL.

### Step 6) Run matrixQTL
This file is a wrapper of file m06_matrixQTL.R.

	python m06_run_matrixQTL.py \
		-p work_path           # same as previous steps

The output of this step is a file named **qtl.txt.gz** which stores the matrixQTL output.

### Step 7) Get cell type specific QTL results
This step split the QTL results for each cell type, because we need to do clump for each cell type separately.

	python m07_cell_qtl.py \
		-p work_path          # same as previous steps
This step creates a folder **cell_qtl** in folder matrixQTL, it splits the QTL results in file qtl.txt.gz into cell type specific QTL results with cell type as the file names.

### Step 8) Q-Q plot checking
	python m08_qqplot.py \
		-p work_path          # same as previous steps
This step generates QQ plot for each cell type QTL results. The output figures are in the same path (cell_qtl) as those QTL result files.

### Step 9) Clump the QTL results for each cell type 
This step clum the QTL for each cell type separately.  
	
	python m09_clump_qtl.py \
		-p work_path \       # same as previous steps
		--pval pvalue        # threshold of pvalue
		--plink2 path/to/plink2/binary \
        --plink1 path/to/plink1/binary		
It generates a folder **clump** in folder matrixQTL. The folder stores clump results for each cell type

### Step 10) Annotate significant QTL results
This step extract significant QTL results and annotate them. You need to have a annotation file generated by VEP. The file needs to have these several columns: 'chr,pos,ref,alt,Consequence,SYMBOL,near_50Kgene,CADD_PHRED,gnomAD_AF,Existing_variation'

	python m10_annotate_clump_qtl.py \
		-p work_path \       # same as previous steps
		--pval pvalue \      # pavalue threshold for significant QTL resutls
		-a anno_fn           # annotated results using VEP


