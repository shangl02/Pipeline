#!/usr/bin/env nextflow

/* ----------- define parameters ---------------------
*/
// (1) related files
ref_fa   =  file(params.genome_fa)
ref_fai  =  ref_fa + '.fai'
ref_dict =  ref_fa.parent / ref_fa.baseName + '.dict'
interval = file(params.interval)

gvcf_path = file(params.gvcf_path)

// (2) check environment resource
ava_mem = (double) (Runtime.getRuntime().freeMemory())
ava_cpu = Runtime.getRuntime().availableProcessors()

if (params.cpu != null && ava_cpu > params.cpu) {
    ava_cpu = params.cpu
}

if (params.mem != null && ava_mem > params.mem) {
    ava_mem = params.mem
}

// (3) reference variation database
dbsnp = file(params.dbsnp)
dbsnp_idx = dbsnp + '.idx'
if (!dbsnp.exists()) exit 1, "dbsnp doesn't exist"

gold_indel = file(params.gold_indel)
gold_indel_idx = gold_indel + '.idx'
if (!gold_indel.exists()) exit 1, "gold indel doesn't exist"

hapmap = file(params.hapmap)
hapmap_idx = hapmap + '.idx'
if (!hapmap.exists()) exit 1, "hapmap vcf doesn't exist"

omni = file(params.omni)
omni_idx = omni + '.idx'
if (!omni.exists()) exit 1, "omni vcf doesn't exist"

phase1 = file(params.phase1)
phase1_idx = phase1 + '.idx'
if (!phase1.exists()) exit 1, "phase1 vcf doesn't exist"

/*----------------- Consolidate variants ------------
*/
Channel
	.fromPath(interval)
	.splitText()
	.map{it.split()[0] + ':' + it.split()[1] + '-' + it.split()[2]}
	.set {inters}


Channel
	.fromPath(interval)
	.splitText()
	.map{it.split()[0] + '_' + it.split()[1] + '_' + it.split()[2]}
	.set {inter_folders}

gvcfs = params.gvcf_path + '/' + params.gvcf_pattern
gvcf_fns = file(gvcfs).collect{"-V $it"}.join(' ')

inters
	.merge(inter_folders)
	.combine(Channel.from(gvcf_fns))
	.set {inter_gvcf}


process ConsolidateGVCF_GenotypeGVCF {
	cache 'deep'

    tag {GVCF_DB}

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    set val(inter), val(inter_folder), val(fn) from inter_gvcf

    output:
    file "${inter_folder}.vcf.gz" into inter_vcf
    file "${inter_folder}.vcf.gz.tbi" into inter_vcf_idx


    script:
    """
    gatk GenomicsDBImport \
    ${fn} \
    --genomicsdb-workspace-path ${inter_folder} \
    --tmp-dir=${gvcf_path} \
    -L ${inter}

    mkdir -p ${params.out_dir}/HaplotypeCaller

    gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R ${ref_fa} \
    -V gendb://${inter_folder} \
    -O ${inter_folder}.vcf.gz \
    --tmp-dir=${params.out_dir}/HaplotypeCaller
    """
}


process mergeGenotypeGVCF {
	publishDir {params.out_dir + '/GenotypeGVCF'}
	
	input:
	file inter_vcf from inter_vcf.collect()
	file inter_vcf_idx from inter_vcf_idx.collect()

	output:
	file "haplotype.vcf.gz" into haplotypecaller_vcf
	file "haplotype.vcf.gz.tbi" into haplotypecaller_vcf_idx
	
	script:
	vcfs = inter_vcf.join(' ')
	"""
	bcftools concat -a -o unsort.vcf.gz -O z ${vcfs}
    mkdir tmp
    bcftools sort -T tmp -o haplotype.vcf.gz -O z unsort.vcf.gz 
    rm unsort.vcf.gz
    tabix haplotype.vcf.gz
    """
}


/*---------------- VQSR --------------------------------
*/
process VQSR_SNP {
    tag {vqsr_snp}
    publishDir pattern: "recal_snp.*",
        path: {params.out_dir + '/VQSR'},
        mode: 'copy', overwrite: true

    input:
    file ref_fa
    file ref_fai
    file ref_dict

    file hapmap
    file hapmap_idx
    file omni
    file omni_idx
    file phase1
    file phase1_idx
    file dbsnp
    file dbsnp_idx

    file haplotypecaller_vcf
    file haplotypecaller_vcf_idx
    file interval

    output:
    file 'recal_snp.recal' into vqsr_snp_recal
    file 'recal_snp.recal.idx' into vqsr_snp_recal_idx
    file 'recal_snp.tranches' into vqsr_snp_tranches
    file 'recal_snp.plots.R' into vqsr_snp_Rscript

    script:
    """
    gatk VariantRecalibrator \
    -R ${ref_fa} \
    -V ${haplotypecaller_vcf} \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${phase1} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
    -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O recal_snp.recal \
    --tranches-file recal_snp.tranches \
    --rscript-file recal_snp.plots.R \
    -L ${interval}
    """
}

process Apply_VQSR_SNP {
    tag {apply_vqsr_snp}

    input:
    file vqsr_snp_recal
    file vqsr_snp_recal_idx
    file vqsr_snp_tranches
    file haplotypecaller_vcf
    file haplotypecaller_vcf_idx
    file interval

    output:
    file 'recal_snps_raw_indels.vcf.gz' into recal_snps_raw_indels
    file 'recal_snps_raw_indels.vcf.gz.tbi' into recal_snps_raw_indels_idx

    script:
    """
    gatk ApplyVQSR \
    -V ${haplotypecaller_vcf} \
    --recal-file ${vqsr_snp_recal} \
    --tranches-file ${vqsr_snp_tranches} \
    -mode SNP \
    -O recal_snps_raw_indels.vcf.gz \
    -L ${interval}
    """

}

process VQSR_INDEL {
    tag {vqsr_indel}
    publishDir pattern: "recal_indel.*",
        path: {params.out_dir + '/VQSR'},
        mode: 'copy', overwrite: true

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    file recal_snps_raw_indels
    file recal_snps_raw_indels_idx

    file dbsnp
    file dbsnp_idx
    file gold_indel
    file gold_indel_idx
    file interval

    output:
    file 'recal_indel.recal' into vqsr_indel_recal
    file 'recal_indel.recal.idx' into vqsr_indel_recal_idx
    file 'recal_indel.tranches' into vqsr_indel_tranches
    file 'recal_indel.plots.R' into vqsr_indel_Rscript
    

    script:
    """
    gatk VariantRecalibrator \
    -V ${recal_snps_raw_indels} \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 ${gold_indel} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
    -an QD -an DP -an FS -an SOR \
    -an MQRankSum -an ReadPosRankSum -mode INDEL \
    --max-gaussians 4 \
    -O recal_indel.recal \
    --tranches-file recal_indel.tranches \
    --rscript-file recal_indel.plots.R \
    -L ${interval}
    """
}

process Apply_VQSR_INDEL {
    tag {apply_vqsr_indel}
    publishDir pattern: "recal_varis.vcf.gz",
        path: {params.out_dir + '/Final'},
        mode: 'copy', overwrite: true

    input:
    file recal_snps_raw_indels
    file recal_snps_raw_indels_idx
    file vqsr_indel_recal
    file vqsr_indel_recal_idx
    file vqsr_indel_tranches
    file interval

    output:
    file 'recal_varis.vcf.gz' into recal_varis_vcf

    script:
    """
    gatk ApplyVQSR \
    -V ${recal_snps_raw_indels} \
    --recal-file ${vqsr_indel_recal} \
    --tranches-file ${vqsr_indel_tranches} \
    -mode INDEL \
    -ts-filter-level 99.0 \
    -O recal_varis.vcf.gz \
    -L ${interval}
    """
}

