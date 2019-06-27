#!/usr/bin/env nextflow

/* ----------- define parameters ---------------------
*/
// (1) related files
ref_fa   =  file(params.genome_fa)
ref_fai  =  ref_fa + '.fai'
ref_dict =  ref_fa.parent / ref_fa.baseName + '.dict'
interval = file(params.interval)

if (!ref_fa.exists()) exit 1, "Reference genome not found: ${ref_fa}"

fq_path = file(params.fq_path)
multiqc_config = file(params.multiqc_config)

// (2) check environment resource
ava_mem = (double) (Runtime.getRuntime().freeMemory())
ava_cpu = Runtime.getRuntime().availableProcessors()

if (params.cpu != null && ava_cpu > params.cpu) {
    ava_cpu = params.cpu
}

if (params.mem != null && ava_mem > params.mem) {
    ava_mem = params.mem
}

idv_cpu = 3
parallel_each_batch = ava_cpu / idv_cpu
if (parallel_each_batch < 1) {
    parallel_each_batch = 1
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
/*----------------- FastQC --------------------------
*/
reads = params.fq_path + '/' + params.fq_pattern

Channel.fromFilePairs(reads, size: params.pair_end ? 2: 1)
        .ifEmpty {
            exit 1, "Can't find fq files"
        }
        .into {reads_for_fastqc; reads_for_trim}

process Run_fastQC {
    tag {fastq_tag}

    publishDir pattern: "*.html",
        path: {params.out_dir + '/Result/QC'},
        mode: 'copy', overwrite: true

    input:
    set val(sample), file(fq_files) from reads_for_fastqc

    output:
    file "*.html" into fastqc_for_waiting

    when:
    params.QC

    script:
    fastq_tag = sample
    fastq_threads = idv_cpu
    fq_file = fq_files.join(' ')
    """
    fastqc -t ${fastq_threads} ${fq_file}
    """
    
}

/*----------------- trim reads --------------------
*/
process Trim_reads {
    tag {trim_tag}

    input:
    set val(sample), file(fq_file) from reads_for_trim
    file tempfiles from fastqc_for_waiting.toList()

    output:
    set val(sample), file(trim_fq_file) into trim_reads_for_bwa, trim_reads_for_index
    file "temp" into trim_for_waiting

    script:
    trim_tag = sample
    trim_fq_file = fq_file.collect {"trim_$it"}
    if (params.trim) {        
        if (params.pair_end) {
            """
            trimmomatic PE -threads ${idv_cpu} -phred33 ${fq_file[0]} ${fq_file[1]} ${trim_fq_file[0]} unpair1.fq.gz ${trim_fq_file[1]} unpair2.fq.gz  ILLUMINACLIP:${params.adapter}:2:30:10 LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
            echo "finish trim" > temp
            """
        } else {
            """
            trimmomatic SE -threads ${idv_cpu} -phred33 ${fq_file} ${trim_fq_file} ILLUMINACLIP:${params.adapter}:2:30:10 LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
            echo "finish trim" > temp
            """
        }
    }
    else {
        if (params.pair_end) {
            """
            mv ${fq_file[0]} ${trim_fq_file[0]}
            mv ${fq_file[1]} ${trim_fq_file[1]}
            echo "finish trim" > temp
            """
        } else {
            """
            mv ${fq_file} ${trim_fq_file}
            echo "finish trim" > temp
            """
        }
    }   
}

/*----------------- align reads using bwa -------------------
*/
// (1) build index
if (params.bwa_index){
    bwa_index = file(params.bwa_index)
    bwa_dir =  bwa_index.parent
    bwa_base = bwa_index.baseName
    if (!bwa_dir.exists()) bwa_dir.mkdir()
    all_idx_files = bwa_dir.list()

}

trim_for_waiting = trim_for_waiting.first()
process Build_index {
    tag {bwa_index}

    publishDir pattern: "*.{ann,amb,bwt,pac,sa}",
        path: {bwa_dir},
        mode: 'copy', overwrite: true

    input:
    file ref_fa
    val bwa_idx from bwa_index
    file tempfiles from trim_for_waiting.toList()

    output:
    file("*.{ann,amb,bwt,pac,sa}") into index_for_waiting

    when:
    all_idx_files.length < 5

    script:
    """
    bwa index -a bwtsw -p ${bwa_idx} ${ref_fa}
    """
}

// (2) align reads using bwa mem
process Align_reads {
    tag {bwa_mem}

    input:
    val bwa_idx from bwa_index
    set val(sample),  file(fq_files) from trim_reads_for_bwa
    file tempfiles from index_for_waiting.toList()

    output:
    set val(sample), file("${sample}_sort.bam") into bam_for_rmdup

    script:
    bwa_threads = ava_cpu
    println fq_files
    if (params.pair_end) {
        """
        bwa mem -M -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina' -t ${bwa_threads} ${bwa_idx} ${fq_files[0]} ${fq_files[1]} | samtools sort -@ ${bwa_threads} -T ${sample} -o ${sample}_sort.bam -
        """
    } else {
        """
        bwa mem -M -R '@RG\\tID:${sample}\\tSM:${sample}\\
        tPL:Illumina' -t ${bwa_threads} ${bwa_idx} ${fq_files} | samtools sort -@ ${bwa_threads} -T ${sample} -o ${sample}_sort.bam -
        """
    }    
}


/*----------------- remove duplicates -------------------
*/
process MarkDuplicates {
    tag {rmdup}

    input:
    set val(sample), file(bam) from bam_for_rmdup

    output:
    set val(sample), file("${sample}_rmdup.bam") into rmdup_bam_for_BQSR_table

    script:
    """
    picard MarkDuplicates I=${bam} M=metrics.txt O=${sample}_rmdup.bam TMP_DIR=tmp
    samtools index ${sample}_rmdup.bam
    """
}

/*---------------- Base recalibration ---------------------
*/
process BuildBQSRTable {
    tag {bqsr_table}

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    file dbsnp
    file dbsnp_idx
    file gold_indel
    file gold_indel_idx
    file interval
    
    set val(sample), file(rmdup_bam) from rmdup_bam_for_BQSR_table

    output:
    set val(sample), file(rmdup_bam) ,file("${sample}.table") into BQSR_table

    script:
    """
    gatk BaseRecalibrator \
    -I ${rmdup_bam} \
    --known-sites ${dbsnp} \
    --known-sites ${gold_indel} \
    -O ${sample}.table \
    -R ${ref_fa} \
    -L ${interval}
    """
}

process RUN_BQSR {
    tag {run_BQSR}

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    set val(sample), file(bam), file(table) from BQSR_table
    file interval

    output:
    set val(sample), file("${sample}_bqsr.bam") into bam_bqsr_for_vari_call

    script:
    """
    gatk ApplyBQSR \
    -R ${ref_fa} \
    -I ${bam} \
    --bqsr-recal-file ${table} \
    -O ${sample}_bqsr.bam \
    -L ${interval}
    samtools index ${sample}_bqsr.bam
    """
}

/*--------------- Call Variants ------------------
*/
process HaplotypeCaller {
    publishDir "${params.out_dir}/HaplotypeCaller"

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    set val(sample), file(bqsr_bam) from bam_bqsr_for_vari_call
    file interval

    output:
    set val(sample), file("${sample}.g.vcf") into vcf_for_genotype
    file "*.g.vcf*" into haplotype_for_waiting

    script:
    """
    gatk HaplotypeCaller -I ${bqsr_bam} \
    -O ${sample}.g.vcf \
    --emit-ref-confidence GVCF \
    -R ${ref_fa} \
    -L ${interval}
    """
}

/*----------------- Consolidate variants ------------
*/
process ConsolidateGVCF {
    tag {mergeGVCF}

    input:
    val gvcf_path from "${params.out_dir}/HaplotypeCaller"
    file tempfiles from haplotype_for_waiting.toList()
    file interval

    output:
    file "genomicsDB" into genomicsDB_path

    script:
    gvcf = file(gvcf_path + '/*.g.vcf').collect{"-V $it"}.join(' ')
    """
    gatk GenomicsDBImport \
    ${gvcf} \
    --genomicsdb-workspace-path genomicsDB \
    --tmp-dir=${gvcf_path} \
    -L ${interval}
    """
}

/*----------------- Run GenotypeGVCFs -----------------------
*/
process GenotypeGVCFs {
    publishDir "${params.out_dir}/GenotypeGVCF"
    
    input:
    file ref_fa
    file ref_fai
    file ref_dict
    file db_path from genomicsDB_path

    output:
    file 'haplotypecaller.vcf' into haplotypecaller_vcf
    file 'haplotypecaller.vcf.idx' into haplotypecaller_vcf_idx
    
    script:
    """
    gatk GenotypeGVCFs \
    -R ${ref_fa} \
    -V gendb://${db_path} \
    -O haplotypecaller.vcf \
    --tmp-dir=${params.out_dir}/HaplotypeCaller
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
    file 'recal_snps_raw_indels.vcf' into recal_snps_raw_indels
    file 'recal_snps_raw_indels.vcf.idx' into recal_snps_raw_indels_idx

    script:
    """
    gatk ApplyVQSR \
    -V ${haplotypecaller_vcf} \
    --recal-file ${vqsr_snp_recal} \
    --tranches-file ${vqsr_snp_tranches} \
    -mode SNP \
    -O recal_snps_raw_indels.vcf \
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
    publishDir pattern: "recal_varis.vcf",
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
    file 'recal_varis.vcf' into recal_varis_vcf

    script:
    """
    gatk ApplyVQSR \
    -V ${recal_snps_raw_indels} \
    --recal-file ${vqsr_indel_recal} \
    --tranches-file ${vqsr_indel_tranches} \
    -mode INDEL \
    -ts-filter-level 99.0 \
    -O recal_varis.vcf \
    -L ${interval}
    """
}
