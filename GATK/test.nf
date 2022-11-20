#!/usr/bin/env nextflow

/* ----------- define parameters ---------------------
*/
// (1) related files
ref_fa   =  file(params.genome_fa)
ref_fai  =  ref_fa + '.fai'
ref_dict =  ref_fa.parent / ref_fa.baseName + '.dict'

if (!ref_fa.exists()) exit 1, "Reference genome not found: ${ref_fa}"

freq_vcf     = file(params.population_freq)
freq_vcf_idx = freq_vcf + '.tbi'

vcfs = params.vcf_path + '/*.vcf.gz'
vcfs_index = params.vcf_path + '/*.vcf.gz.tbi'
peds = params.vcf_path + '/*.ped'

vcf_files = Channel.fromPath(vcfs, checkIfExists: true)
				.map {file -> [file.simpleName, file]}

vcf_indexes = Channel.fromPath(vcfs_index, checkIfExists: true)
				.map {file -> [file.simpleName, file]}

ped_files = Channel.fromPath(peds, checkIfExists: true)
				.map {file -> [file.simpleName, file]}


vcf_files
	.combine(vcf_indexes, by: 0)
	.combine(ped_files, by: 0)
	.set {vcf_for_post}   // format: [family, vcf, vcf_idx, ped]

/* ------------- CalculateGenotypePosteriors -------
*/
process CalculateGenotypePosteriors {
	tag "${family}"
	
	input:
	set val(family), file(vcf), file(vcf_idx), file(ped) from vcf_for_post
	file(freq_vcf)
	file(freq_vcf_idx)

	output:
	set val(family), file("${family}.GT_post.vcf.gz"), file("${family}.GT_post.vcf.gz.tbi"), file(ped) into gt_post_for_filter, gt_for_merge

	script:
	"""	
	gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
	-V ${vcf} -O ${family}.GT_post.vcf.gz \
	-ped ${ped} \
	-supporting ${freq_vcf}

	echo ${family}.vcf.gz.tbi
	echo ${ped}
	"""
}

Channel
    .from(1..22)
    .map {'chr' + it.toString()}
    .set {chr_idx}

other = Channel.from('chrX','chrY','chrM')
chr_idx.concat(other).set {chroms}

chroms
    .combine(gt_post_for_filter)
    .set {chrom_vcf}  // format: [chr, family, vcf, vcf_idx, ped]


/* ----------- Filter low GQ ----------------------
*/
process VariantFilter {
	tag "${family}.${chr}"

	input:
	set val(chr), val(family), file(vcf), file(vcf_idx), file(ped) from chrom_vcf
	file(ref_fa)
	file(ref_fai)
	file(ref_dict)

	output:
	file("${chr}.rm_lowGQ.vcf.gz") into chrom_g_vcf
	file("${chr}.rm_lowGQ.vcf.gz.tbi") into chrom_g_vcf_index
	file "*.vcf.gz" into waiting

	script:
	"""
	gatk VariantFiltration \
	   -R ${ref_fa} -V ${vcf} \
	   -O ${chr}.lowGQ.vcf.gz \
	   -G-filter-name "lowGQ" \
	   -G-filter "GQ < 20" \
	   -L ${chr}

	gatk SelectVariants \
		-V ${chr}.lowGQ.vcf.gz \
		--set-filtered-gt-to-nocall \
		-O ${chr}.rm_lowGQ.vcf.gz

	echo ${chr}.rm_lowGQ.vcf.gz.tbi
	"""
}


process merge {
	
	input:
	file chroms from chrom_g_vcf.collect()
	file index from chrom_g_vcf_index.collect()
	set val(family), file(vcf), file(idx), file(ped) from gt_for_merge

	script:
	"""
	echo test
	"""


}