#!/usr/bin/env nextflow

/* ----------- define parameters ---------------------
*/
// (1) check environment resource
ava_mem = (double) (Runtime.getRuntime().freeMemory())
ava_cpu = Runtime.getRuntime().availableProcessors()

if (params.cpu != null && ava_cpu > params.cpu) {
    ava_cpu = params.cpu
}

if (params.mem != null && ava_mem > params.mem) {
    ava_mem = params.mem
}


// (2) related files
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
	.set {vcf_and_ped}   // format: [sample, vcf, vcf_idx, ped]

// get all of the chromosomes
Channel
    .from(1..22)
    .map {'chr' + it.toString()}
    .set {chr_idx}

other = Channel.from('chrX','chrY','chrM')
chr_idx
	.concat(other)
	.into {chroms; chroms_for_merge}

// cross all vcf files with all chromosomes
vcf_and_ped
	.combine(chroms)
	.set {vcf_for_post} // format: [sample,vcf,vcf_idx,ped,chrom]


/* ------------- CalculateGenotypePosteriors -------
*/
process CalculateGenotypePosteriors {
	tag "${family}.${chr}"
	
	input:
	set val(family), file(vcf), file(vcf_idx), file(ped), \
	val(chr) from vcf_for_post
	file(freq_vcf)
	file(freq_vcf_idx)

	output:
	set val(family), val(chr), \
	file("${family}.${chr}.GT_post.vcf.gz"), \
	file("${family}.${chr}.GT_post.vcf.gz.tbi"), \
	file(ped) into gt_post_for_filter

	script:
	"""	
	gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
	-V ${vcf} -O ${family}.${chr}.GT_post.vcf.gz \
	-ped ${ped} -L ${chr} \
	-supporting ${freq_vcf}

	echo ${family}.${chr}.vcf.gz.tbi
	echo ${ped}
	"""
}


/* ----------- Filter low GQ ----------------------
*/
process VariantFilter {
	tag "${family}.${chr}"

	input:
	set val(family), val(chr), file(vcf), file(vcf_idx), \
	file(ped) from gt_post_for_filter
	file(ref_fa)
	file(ref_fai)
	file(ref_dict)

	output:
	set val(family), val(chr), \
	file("${family}.${chr}.rm_lowGQ.vcf.gz"), \
	file("${family}.${chr}.rm_lowGQ.vcf.gz.tbi"), \
	file(ped) into rm_lowGQ_for_denovo

	script:
	"""
	gatk VariantFiltration \
	   -R ${ref_fa} -V ${vcf} \
	   -O ${family}.${chr}.lowGQ.vcf.gz \
	   -G-filter-name "lowGQ" \
	   -G-filter "GQ < 20"

	gatk SelectVariants \
		-V ${family}.${chr}.lowGQ.vcf.gz \
		--set-filtered-gt-to-nocall \
		-O ${family}.${chr}.rm_lowGQ.vcf.gz

	echo ${family}.${chr}.rm_lowGQ.vcf.gz.tbi
	"""
}


/* -------------- find denovo -------------
*/
process Denovo {
	tag "${family}.${chr}"

	publishDir pattern: "*.{gz,tbi}",
        path: {params.out_dir + "/denovo/${family}"},
        mode: 'copy', overwrite: true

	input:
	set val(family), val(chr), file(vcf), file(vcf_idx), \
	file(ped) from rm_lowGQ_for_denovo
	file(ref_fa)
	file(ref_fai)
	file(ref_dict)

	output:
	set val(family), val(chr), file("${family}.${chr}.denovo.vcf.gz"), file("${family}.${chr}.denovo.vcf.gz.tbi") into denovo_vcf

	script:
	"""
	gatk VariantAnnotator \
    -R ${ref_fa} -V ${vcf} \
    -O ${family}.${chr}.denovo.vcf.gz \
    -A PossibleDeNovo \
    -ped ${ped}
    echo ${family}.${chr}.denovo.vcf.gz.tbi
	"""
}


/* ------------------ Merge vcf files for families ----------
*/

denovo_vcf   // format: [family, chr, vcf, vcf_idx]
	.groupTuple()
	.set {all_vcfs}

process MergeVCF {
	tag "${family}"

	publishDir pattern: "*.{gz,tbi}",
        path: {params.out_dir + "/denovo"},
        mode: 'copy', overwrite: true

    input:
    set val(family), val(chroms), file(vcfs), file(vcf_idx) from all_vcfs

    output:
    set val(family), file("${family}.denovo.vcf.gz"), \
    file("${family}.denovo.vcf.gz.tbi") into final_denovo

    script:
    """
    bcftools concat -o ${family}.denovo.vcf.gz -O z --threads ${ava_cpu} ${family}.chr1.denovo.vcf.gz \
    ${family}.chr2.denovo.vcf.gz \
    ${family}.chr3.denovo.vcf.gz \
    ${family}.chr4.denovo.vcf.gz \
    ${family}.chr5.denovo.vcf.gz \
    ${family}.chr6.denovo.vcf.gz \
    ${family}.chr7.denovo.vcf.gz \
    ${family}.chr8.denovo.vcf.gz \
    ${family}.chr9.denovo.vcf.gz \
    ${family}.chr10.denovo.vcf.gz \
    ${family}.chr11.denovo.vcf.gz \
    ${family}.chr12.denovo.vcf.gz \
    ${family}.chr13.denovo.vcf.gz \
    ${family}.chr14.denovo.vcf.gz \
    ${family}.chr15.denovo.vcf.gz \
    ${family}.chr16.denovo.vcf.gz \
    ${family}.chr17.denovo.vcf.gz \
    ${family}.chr18.denovo.vcf.gz \
    ${family}.chr19.denovo.vcf.gz \
    ${family}.chr20.denovo.vcf.gz \
    ${family}.chr21.denovo.vcf.gz \
    ${family}.chr22.denovo.vcf.gz \
    ${family}.chrX.denovo.vcf.gz \
    ${family}.chrY.denovo.vcf.gz \
    ${family}.chrM.denovo.vcf.gz

    tabix ${family}.denovo.vcf.gz
    echo ${family}.denovo.vcf.gz.tbi
    """
}





