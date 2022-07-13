### This pipeline runs finemapping using SusieR.
You can download the docker container in HPC using the following command:
	
	singularity pull docker://shl198/finemap_susie:0.11.97

#### Some points need to pay attention
* input file needs to be tab separated.
* for --stats parameter, if the column name has '-' in it, change it to '.' in the parameter, eg: columns name is 't-stats', then the prameter is --stats t.stats.


#### How to run the pipeline

	singularity run -B /:/media finemap.sif \
		python3 /opt/finemapping.py \
		-i /path/to/gwas/file \
		-r /path/to/reference/genotype/file/in/plink/format \
		-s SNP \
		-o /output/path \
		--stats z_score

* -i: input gwas file, can be gzip or not, needs to be full path.
* -r: reference genotype file for calculating LD in plink format, needs to be the binary file which ends with bed/bgen.
* -s: column name indicating SNP id in the input gwas file.
* -o: output path that will save all the results
* --stats: column name in gwas file indicating the stats of the results

