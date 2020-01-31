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

rg_file = file(params.RG_file)
if (!rg_file.exists()) exit 1, "ReadGroup design file not found: ${rg_file}"

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
if (params.pair_end) {
    Channel
        .fromPath(rg_file)
        .splitCsv(header:true, sep:'\t')
        .map {row -> [row.sample_id, row.lane_id, [file(params.fq_path+'/'+row.fq_1), file(params.fq_path+'/'+row.fq_2)]]}
        .into {reads_for_fastqc; reads_for_trim}
} else {
    Channel
        .fromPath(rg_file)
        .splitCsv(header:true, sep:'\t')
        .map {row -> [row.sample_id, row.lane_id, [file(params.fq_path+'/'+row.fq_1)]]}
        .into {reads_for_fastqc; reads_for_trim}

}


process Run_fastQC {
    tag "${sample}_${lane}"

    publishDir pattern: "*.{html,zip}",
        path: {params.out_dir + '/MultiQC/QC01_fastqc'},
        mode: 'copy', overwrite: true

    input:
    set val(sample), val(lane), file(fq_files) from reads_for_fastqc

    output:
    file "*.html" into fastqc_for_waiting
    file "*.zip" into fastqc_zip

    when:
    params.QC

    script:
    if (params.pair_end) {
        """
        [ ! -f ${sample}_${lane}_1.fq.gz ] && ln -s ${fq_files[0]} ${sample}_${lane}_1.fq.gz
        [ ! -f ${sample}_${lane}_2.fq.gz ] && ln -s ${fq_files[1]} ${sample}_${lane}_2.fq.gz
        fastqc -t 2 ${sample}_${lane}_1.fq.gz ${sample}_${lane}_2.fq.gz
        """
    } else {
        """
        [ ! -f ${sample}_${lane}.fq.gz ] && ln -s ${fq_files[0]} ${sample}_${lane}.fq.gz
        fastqc ${sample}_${lane}.fq.gz
        """
    }
}


/*----------------- trim reads --------------------
*/
if (! params.trim) {
    reads_for_trim.into {
        trim_reads_for_bwa;
        trim_reads_for_index;
        trim_log
    }    
} else {
    process Trim_reads {
        tag "${sample}_${lane}"

        publishDir pattern: "*.log",
            path: {params.out_dir + '/MultiQC/QC02_trimmomatic'},
            mode: 'copy', overwrite: true

        input:
        set val(sample), val(lane), file(fq_file) from reads_for_trim
        file tempfiles from fastqc_for_waiting.toList()


        output:
        set val(sample), val(lane), file(trim_fq_file) into trim_reads_for_bwa, trim_reads_for_index
        file "${sample}_${lane}.log" into trim_log

        script:
        // create link of fq files to make consistant naming in multiQC
        trim_tag = sample
        if (params.pair_end) {
            new_file_1 = "${sample}_${lane}_1.fq.gz"
            new_file_2 = "${sample}_${lane}_2.fq.gz"
            trim_fq_file = ["trim_${new_file_1}", "trim_${new_file_2}"]
            if (params.adapter == '') {
                """
                [ ! -f ${new_file_1} ] && ln -s ${fq_file[0]} ${new_file_1}
                [ ! -f ${new_file_2} ] && ln -s ${fq_file[1]} ${new_file_2}
                trimmomatic PE -threads ${ava_cpu} -phred33 ${new_file_1} ${new_file_2} trim_${new_file_1} unpair1.fq.gz trim_${new_file_2} unpair2.fq.gz LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 2> ${sample}_${lane}.log
                """
            } else {
                """
                [ ! -f ${new_file_1} ] && ln -s ${fq_file[0]} ${new_file_1}
                [ ! -f ${new_file_2} ] && ln -s ${fq_file[1]} ${new_file_2}
                trimmomatic PE -threads ${ava_cpu} -phred33 ${new_file_1} ${new_file_2} trim_${new_file_1} unpair1.fq.gz trim_${new_file_2} unpair2.fq.gz ILLUMINACLIP:${params.adapter}:2:30:10 LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 2> ${sample}_${lane}.log
                """
            }
        } else { // single_end
            new_file = "${sample}_${lane}.fq.gz"
            trim_fq_file = ["trim_${new_file}"]
            if (params.adapter == '') {
                """
                [ ! -f ${new_file} ] && ln -s ${fq_file[0]} ${new_file}
                trimmomatic SE -threads ${ava_cpu} -phred33 ${new_file} trim_${new_file} LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 2> ${sample}_${lane}.log
                """
            } else {
                """
                [ ! -f ${new_file} ] && ln -s ${fq_file[0]} ${new_file}
                trimmomatic SE -threads ${ava_cpu} -phred33 ${new_file} trim_${new_file} ILLUMINACLIP:${params.adapter}:2:30:10 LEADING:15 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 2> ${sample}_${lane}.log
                """
            }
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

trim_for_waiting = trim_log.first()
process Build_index {

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
    tag "${sample}_${lane}"

    publishDir pattern: "*.log",
        path: {params.out_dir + '/MultiQC/QC03_bwaBam'},
        mode: 'copy', overwrite: true

    input:
    val bwa_idx from bwa_index
    set val(sample), val(lane), file(fq_files) from trim_reads_for_bwa
    file tempfiles from index_for_waiting.toList()

    output:
    set val(sample), file("${sample}_${lane}_sort.bam") into bam_for_merge
    set val(sample), file("${sample}_${lane}_sort.bam.bai") into bam_for_merge_index
    file "${sample}_${lane}.log" into bam_stats


    script:
    bwa_threads = ava_cpu
    if (params.pair_end) {
        """
        bwa mem -M -R '@RG\\tID:${sample}_${lane}\\tSM:${sample}\\tPL:Illumina' -t ${bwa_threads} ${bwa_idx} ${fq_files[0]} ${fq_files[1]} | samtools sort -@ ${bwa_threads} -T ${sample}_${lane} -o ${sample}_${lane}_sort.bam -
        samtools index ${sample}_${lane}_sort.bam
        samtools flagstat ${sample}_${lane}_sort.bam > ${sample}_${lane}.log
        """
    } else {
        """
        bwa mem -M -R '@RG\\tID:${sample}_${lane}\\tSM:${sample}\\        tPL:Illumina' -t ${bwa_threads} ${bwa_idx} ${fq_files} | samtools sort -@ ${bwa_threads} -T ${sample}_${lane} -o ${sample}_${lane}_sort.bam -
        samtools index ${sample}_${lane}_sort.bam
        samtools flagstat ${sample}_${lane}_sort.bam > ${sample}_${lane}.log
        """
    }    
}

bam_for_merge
    .groupTuple()
    .set {merged_bam}

bam_for_merge_index
    .groupTuple()
    .set {merged_bam_index}

/*----------------- Merge Bams from different lanes ---------
*/
process MergeBam {
    tag "${sample}"

    publishDir pattern: "*.log",
        path: {params.out_dir + '/MultiQC/QC04_mergeBam'},
        mode: 'copy', overwrite: true

    input:
    set val(sample), file(bams) from merged_bam
    set val(sample), file(bams_index) from merged_bam_index

    output:
    set val(sample), file("${sample}.sort.bam") into bam_for_rmdup
    file("${sample}.log") into merged_bam_log

    script:
    merge_threads = ava_cpu
    bam_files = bams.join(' ')
    """
    samtools merge -@ ${merge_threads} ${sample}.sort.bam ${bam_files}
    samtools index ${sample}.sort.bam
    samtools flagstat ${sample}.sort.bam > ${sample}.log
    """
}

/*----------------- remove duplicates -------------------
*/
process MarkDuplicates {
    tag "${sample}"

    publishDir pattern: "*.log",
        path: {params.out_dir + '/MultiQC/QC05_MarkDuplicates'},
        mode: 'copy', overwrite: true

    input:
    set val(sample), file(bam) from bam_for_rmdup

    output:
    set val(sample), file("${sample}_rmdup.bam") into rmdup_bam_for_BQSR_table
    file "${sample}.log" into markdup_log

    script:
    """
    gatk --java-options "-Xmx8G" MarkDuplicates --INPUT ${bam} --METRICS_FILE ${sample}.log --OUTPUT ${sample}_rmdup.bam --TMP_DIR tmp
    samtools index ${sample}_rmdup.bam
    """
}

/*---------------- Base recalibration ---------------------
*/
process BuildBQSRTable {
    tag "${sample}"

    publishDir pattern: "*.table",
        path: {params.out_dir + '/MultiQC/QC06_BQSR'},
        mode: 'copy', overwrite: true

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
    tag "${sample}"

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    set val(sample), file(bam), file(table) from BQSR_table
    file interval

    output:
    set val(sample), file("${sample}_bqsr.bam") into bam_bqsr_for_vari_call, bam_for_merge_gvcf

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
    .set {chrom_bams}

process HaplotypeCaller {
    tag "${sample}_${chrom}"

    input:
    file ref_fa
    file ref_fai
    file ref_dict
    set chrom, val(sample), file(bqsr_bam) from chrom_bams

    output:
    file "${chrom}.g.vcf.gz" into chrom_g_vcf
    file "${chrom}.g.vcf.gz.tbi" into chrom_g_vcf_index
    file "*.g.vcf*" into haplotype_for_waiting

    script:
    """
    gatk HaplotypeCaller -I ${bqsr_bam} \
    -O ${chrom}.g.vcf \
    --emit-ref-confidence GVCF \
    -R ${ref_fa} \
    -L ${chrom}

    bgzip ${chrom}.g.vcf
    tabix ${chrom}.g.vcf.gz
    echo ${chrom}.g.vcf.gz.tbi
    """
}

/*--------------  Merge variants for all chromosomes --------------------
*/
process MergeVCF {
    tag "${sample}"

    publishDir pattern: "*.g.vcf.gz",
        path: {params.out_dir + '/HaplotypeCaller'},
        mode: 'copy', overwrite: true

    input:
    file chroms from chrom_g_vcf.collect()
    file index from chrom_g_vcf_index.collect()
    set val(sample), file(bam) from bam_for_merge_gvcf

    output:
    file "${sample}.g.vcf.gz" into merge_g_vcf

    script:
    """
    bcftools concat -o ${sample}.g.vcf.gz -O z chrM.g.vcf.gz chr1.g.vcf.gz chr2.g.vcf.gz chr3.g.vcf.gz chr4.g.vcf.gz chr5.g.vcf.gz chr6.g.vcf.gz chr7.g.vcf.gz chr8.g.vcf.gz chr9.g.vcf.gz chr10.g.vcf.gz chr11.g.vcf.gz chr12.g.vcf.gz chr13.g.vcf.gz chr14.g.vcf.gz chr15.g.vcf.gz chr16.g.vcf.gz chr17.g.vcf.gz chr18.g.vcf.gz chr19.g.vcf.gz chr20.g.vcf.gz chr21.g.vcf.gz chr22.g.vcf.gz chrX.g.vcf.gz chrY.g.vcf.gz 
    tabix ${sample}.g.vcf.gz
    """
}
