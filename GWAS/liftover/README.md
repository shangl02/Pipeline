### This docker file has the pipeline to liftover GWAS summary stats
You can download the docker container in HPC using the following command:

	singularity pull docker://shl198/crossmap4gwas:0.6.v2

#### If your GWAS summary stats file has one position column, use the following code:
	
	singularity run -B /:/media crossmap4gwas_0.6.v2.sif \
		python /opt/crossmap_gwas_summary.py \
		-i /lustre/scratch/lis262/crossmap/test37.gwas.gz \
		-o /lustre/scratch/lis262/crossmap/test38.gwas.gz \
		-t /lustre/scratch/lis262/crossmap/temp \
		-c '#CHROM' \
		-p 'POS' \
		--conversion 37to38 \
		--add_end yes

* -i: input gwas file, can be gzip or not, needs to be full path
* -o: output gwas file, output file will be gzip to save space, needs to be full path.
* -t: temparay folder to store some intermediate files, needs to be full path
* -c: header in the input gwas file that indicates chromosome
* -p: header in the input gwas file that indicates position
* --conversion: 38to37 or 37to38
* --add_end: to inidicate if add an end column in the result, the end position will be pos+1. if you don't want that columns, set it to 'no'.

#### If your GWAS summary stats file has chromosome and position columns merged together, use the following code
	singularity run -B /:/media crossmap4gwas_0.6.v2.sif \
		python /opt/crossmap_gwas_summary.py \
		-i /lustre/scratch/lis262/crossmap/test37.gwas.gz \
		-o /lustre/scratch/lis262/crossmap/test38.gwas.gz \
		-t /lustre/scratch/lis262/crossmap/temp \
		-m MarkName \
		--sep ':' \
		--conversion 37to38


* -i: input gwas file, can be gzip or not, needs to be full path
* -o: output gwas file, output file will be gzip to save space, needs to be full path.
* -t: temparay folder to store some intermediate files, needs to be full path
* -m: header in the input gwas file that indicates merged chromosome and position
* --sep: separator for chromosome and position in the marker column
* --conversion: 38to37 or 37to38


#### If your GWAS summary stats file has both start and end column, use the following code:
	
	singularity run -B /:/media crossmap4gwas_0.6.v2.sif \
		python /opt/crossmap_gwas_summary.py \
		-i /lustre/scratch/lis262/crossmap/test37.gwas.gz \
		-o /lustre/scratch/lis262/crossmap/test38.gwas.gz \
		-t /lustre/scratch/lis262/crossmap/temp \
		-c '#CHROM' \
		-s 'pos_start' \
		-e 'pos_end'   \
		--conversion 37to38

* -i: input gwas file, can be gzip or not, needs to be full path
* -o: output gwas file, output file will be gzip to save space, needs to be full path.
* -t: temparay folder to store some intermediate files, needs to be full path
* -c: header in the input gwas file that indicates chromosome
* -s: header in the input gwas file that indicates start position
* -e: header in the input gwas file that indicates end position
* --conversion: 38to37 or 37to38
