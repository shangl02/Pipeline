This pipeline uses sra toolkit fasterq-dump to download raw fastq files.
First, it read in gse desing file, then get srr id from srx id, and use bsub to download each srr fastq files in parallel, after the download it changes file name to gseid.fastq.gz

## Environment preparation

	conda env create -f environment.yml


## ncbi parameter configuration
By default, fasterq-dump download the cache sra file into your home folder, if you are using servers like HPC, you may be allocated limited space, so it is recommended to change the default cache file location. To do that, you need to modify the parameters in file **~/.ncbi/user-settings.mkfg**

	/repository/user/default-path = "/other/path"
	/repository/user/main/public/root = "/other/path"

## Run the pipeline

	python wrapper_download_srr.py -i gene_design_file.xls -o output_path