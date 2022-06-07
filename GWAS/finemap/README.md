### This pipeline runs finemapping using SusieR.

#### Some points need to be 
* input file needs to be tab separated.
* for --stats parameter, if the column name has '-' in it, change it to '.'


#### How to run the pipeline

	singularity run -B /:/media finemap.sif \
		python3 /opt/finemapping.py \
		-i /path/to/gwas/file \
		-r /path/to/reference/genotype/file/in/plink/format \
		-s SNP \
		-o /output/path \
		--stats z_score

* -i: input gwas file, can be gzip or not, needs to be full path. The file needs to be 
* -r: reference genotype file for calculating LD in plink format, needs to be the binary file which ends with bed/bgen.
* -s: column name indicating SNP id in the input gwas file.
* -o: output path that will save all the results
* --stats: column name in gwas file indicating the stats of the results

