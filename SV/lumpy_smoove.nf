#!/usr/bin/env nextflow

/* ----------- define parameters ---------------------
*/
// (1) related files
ref_fa   =  file(params.genome_fa)
ref_fai  =  ref_fa + '.fai'
ref_dict =  ref_fa.parent / ref_fa.baseName + '.dict'

exclude_interval = file(params.exclude_interval)

// (2) check environment resource
ava_mem = (double) (Runtime.getRuntime().freeMemory())
ava_cpu = Runtime.getRuntime().availableProcessors()

if (params.cpu != null && ava_cpu > params.cpu) {
    ava_cpu = params.cpu
}

if (params.mem != null && ava_mem > params.mem) {
    ava_mem = params.mem
}


/*-------------- Dedup bam files -----------------------
*/
bam_files = params.bam_path + '/' + params.bam_pattern
Channel
	.fromPath(bam_files)
	.set{bam_for_dedup}


process BamDedup {
	
	input:
	file bam from bam_for_dedup

	output:
	set val("${bam.baseName}"), file("${bam.baseName}.dedup.bam"), file("${bam.baseName}.dedup.bam.bai") into dedup_bam,bam_for_genotype
	file "${bam.baseName}.metrics.txt" into metrics_for_waiting

	script:
	"""
	picard MarkDuplicates I=${bam} O=${bam.baseName}.dedup.bam M=${bam.baseName}.metrics.txt TMP_DIR=tmp
	samtools index ${bam.baseName}.dedup.bam
	echo ${bam.baseName}.dedup.bam.bai
	"""
}

/*-------------- Run Smoove -----------------------
*/

// step 1: call genotypes

process CallSmoove {
	publishDir pattern: "*.vcf.gz*",
        path: {params.out_dir + '/f01_lumpy_result'},
        mode: 'copy', overwrite: true

	input:
	file ref_fa
	file ref_fai
	file exclude_interval
	set val(sample), file(dedup), file(dedup_idx) from dedup_bam
	file metrics from metrics_for_waiting
	
	output:
	file "${sample}/${sample}-smoove.genotyped.vcf.gz" into sv_vcf
	file "${sample}/${sample}-smoove.genotyped.vcf.gz.csi" into sv_vcf_idx

	script:
	"""
	smoove call --outdir ${sample} --exclude ${exclude_interval} --name ${sample} --fasta ${ref_fa} -p 1 --genotype ${dedup}
	"""
}


// step2: get union set of SV

process MergeSV {
	publishDir pattern: "merged.sites.vcf.gz",
	path: {params.out_dir + '/f02_merged_SV'},
	mode: 'copy', overwrite: true
	
	input:
	file ref_fa
	file ref_fai
	file sv from sv_vcf.collect()
	file sv_idx from sv_vcf_idx.collect()

	output:
	file "merged.sites.vcf.gz" into merged_vcf

	script:
	"""
	smoove merge --name merged -f ${ref_fa} *.genotyped.vcf.gz
	echo merged.sites.vcf.gz
	"""
}

// step3: genotype all samples at the union sites

bam_for_genotype
	.combine(merged_vcf)
	.set {genotype_bam}

process Genotype {
	publishDir pattern: "genotyped/*",
	path: {params.out_dir + '/f03_genotype_SV'},
	mode: 'copy', overwrite: true

	input:
	file ref_fa
	file ref_fai
	set val(sample), file(bam), file(bam_idx), file(vcf) from genotype_bam

	output:
	file "genotyped/${sample}-smoove.genotyped.vcf.gz" into genotyped
	file "genotyped/${sample}-smoove.genotyped.vcf.gz.csi" into genotyped_idx

	script:
	"""
	smoove genotype -d -x -p 1 --name ${sample} --outdir genotyped --fasta ${ref_fa} --vcf ${vcf} ${bam}
	echo genotyped/${sample}-smoove.genotyped.vcf.gz
	echo genotyped/${sample}-smoove.genotyped.vcf.gz.csi
	"""
}

// step4: joint call

process Joint_call {
	publishDir pattern: "*",
	path: {params.out_dir + '/f04_Final_SV'},
	mode: 'copy', overwrite: true

	input:
	file genotype_vcf from genotyped.collect()
	file genotype_vcf_idx from genotyped_idx.collect()

	output:
	file "cohort.smoove.square.vcf.gz"
	file "cohort.smoove.square.vcf.gz.tbi"

	script:
	"""
	smoove paste --name cohort *.vcf.gz
	tabix cohort.smoove.square.vcf.gz
	echo cohort.smoove.square.vcf.gz.tbi
	"""
}