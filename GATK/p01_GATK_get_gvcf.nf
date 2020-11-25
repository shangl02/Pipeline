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

Channel.fromFilePairs(reads, size: params.pair_end ? 2 : 1)
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
    fastq_threads = params.pair_end ? 2 : 1
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
            trimmomatic PE -threads ${ava_cpu} -phred33 ${fq_file[0]} ${fq_file[1]} ${trim_fq_file[0]} unpair1.fq.gz ${trim_fq_file[1]} unpair2.fq.gz  ILLUMINACLIP:${params.adapter}:2:30:10 LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
            echo "finish trim" > temp
            """
        } else {
            """
            trimmomatic SE -threads ${ava_cpu} -phred33 ${fq_file} ${trim_fq_file} ILLUMINACLIP:${params.adapter}:2:30:10 LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
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
    tag {bwa_mem_tag}

    input:
    val bwa_idx from bwa_index
    set val(sample),  file(fq_files) from trim_reads_for_bwa
    file tempfiles from index_for_waiting.toList()

    output:
    set val(sample), file("${sample}_sort.bam") into bam_for_rmdup

    script:
    bwa_mem_tag = sample
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
    rmdup = sample
    """
    gatk --java-options "-Xmx8G" MarkDuplicates --INPUT ${bam} --METRICS_FILE metrics.txt --OUTPUT ${sample}_rmdup.bam --TMP_DIR tmp
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
    bqsr_table = sample
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
    set val(sample), file("${sample}_bqsr.bam") into bam_bqsr_for_vari_call, bam_for_merge_gvcf

    script:
    run_BQSR = sample
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

/*--------------- Call Variants for each chromosome------------------
*/
Channel
    .from(1..22)
    .map {'chr' + it.toString()}
    .set {chr_idx}

other = Channel.from('chrX','chrY','chrM')
chr_idx.concat(other).set {chroms}

chroms
    .combine(bam_bqsr_for_vari_call)
    .set {chrom_bams}    // format: [chrom, sample, bam_file]

process HaplotypeCaller {
    tag "${sample}${chrom}"

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    set val(chrom), val(sample), file(bqsr_bam) from chrom_bams

    output:
    set val(sample), val(chrom), file("${sample}.${chrom}.g.vcf.gz"), file("${sample}.${chrom}.g.vcf.gz.tbi") into chrom_g_vcf
    

    script:
    """
    gatk HaplotypeCaller -I ${bqsr_bam} \
    -O ${sample}.${chrom}.g.vcf.gz \
    --emit-ref-confidence GVCF \
    -R ${ref_fa} \
    -L ${chrom}

    echo ${sample}.${chrom}.g.vcf.gz.tbi
    """
}

/*--------------  Merge variants for all chromosomes --------------------
*/
chrom_g_vcf  // [sample, chr, vcf, vcf_index]
    .groupTuple()
    .set {all_vcfs}


process MergeVCF {
    tag "${sample}"

    publishDir pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}",
        path: {params.out_dir + '/HaplotypeCaller'},
        mode: 'copy', overwrite: true

    input:
    set val(sample), val(chrom), file(vcf), file(vcf_idx) from all_vcfs

    output:
    set file("${sample}.g.vcf.gz"), file("${sample}.g.vcf.gz.tbi") into merge_g_vcf

    script:
    """
    bcftools concat -o ${sample}.g.vcf.gz -O z \
    ${sample}.chr1.g.vcf.gz \
    ${sample}.chr2.g.vcf.gz \
    ${sample}.chr3.g.vcf.gz \
    ${sample}.chr4.g.vcf.gz \
    ${sample}.chr5.g.vcf.gz \
    ${sample}.chr6.g.vcf.gz \
    ${sample}.chr7.g.vcf.gz \
    ${sample}.chr8.g.vcf.gz \
    ${sample}.chr9.g.vcf.gz \
    ${sample}.chr10.g.vcf.gz \
    ${sample}.chr11.g.vcf.gz \
    ${sample}.chr12.g.vcf.gz \
    ${sample}.chr13.g.vcf.gz \
    ${sample}.chr14.g.vcf.gz \
    ${sample}.chr15.g.vcf.gz \
    ${sample}.chr16.g.vcf.gz \
    ${sample}.chr17.g.vcf.gz \
    ${sample}.chr18.g.vcf.gz \
    ${sample}.chr19.g.vcf.gz \
    ${sample}.chr20.g.vcf.gz \
    ${sample}.chr21.g.vcf.gz \
    ${sample}.chr22.g.vcf.gz \
    ${sample}.chrX.g.vcf.gz \
    ${sample}.chrY.g.vcf.gz \
    ${sample}.chrM.g.vcf.gz
    tabix ${sample}.g.vcf.gz
    echo ${sample}.g.vcf.gz.tbi
    """
}
