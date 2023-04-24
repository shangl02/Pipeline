#!/usr/bin/env nextflow

/* This pipeline runs QTL analysis between SNPs and cell type abundance.
*/

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


process QC {
    tag "${plink_pre}"
    publishDir "${params.work_path}/f01_QC", mode: "copy"

    input:
        path vcf_fn
        each path(id_map_fn)

    output:
        tuple val(plink_pre), path("${plink_pre}*"), emit: qc_out

    script:
    plink_pre = vcf_fn.simpleName
        """
        python3 /opt/utils/m01_process_gt.py --vcf ${vcf_fn} --af 0.01 --pre ${plink_pre}
        python3 /opt/utils/m02_clean_samples.py --pre ${plink_pre} --id_map ${id_map_fn}
        """
}

process PCA {
    tag "${plink_pre}"
    publishDir "${params.work_path}/f02_PCA", mode:"copy"

    input:
        tuple val(plink_pre), path(pre_files)
        tuple val(ld_ref_pre), path(ld_ref)
        each path(high_ld)

    output:
        tuple val("${plink_pre}"), path("${plink_pre}*"), emit: pca_out
    
    script:
        """
        python3 /opt/utils/m03_calculate_PC.py --ref ${ld_ref_pre} --hld ${high_ld} --pre ${plink_pre}
        """
}


process get_cova_fn {
    tag "${plink_pre}"
    publishDir "${params.work_path}/f03_Cova", mode:"copy"

    input:
        tuple val(plink_pre), path(pre_files)
        each path(id_map_fn)
        each path(cyto_fn)

    output:
        tuple val(plink_pre), path("${plink_pre}*"), emit: cova_out

    script:
        """
        python3 /opt/utils/m04_get_covariate_file.py --id_map ${id_map_fn} -d ${cyto_fn} --pre ${plink_pre}
        """
}


process matrix_qtl {
    tag "${plink_pre}"
    publishDir "${params.work_path}/f04_matrixQTL", mode:"copy"

    input:
        tuple val(plink_pre), path(plink_files)
        tuple val(pca_pre), path(pca_files)
        tuple val(cova_pre), path(cova_files)

    output:
        tuple val(pca_pre), path("${pca_pre}*"), emit: qtl_out

    script:
        """
        python3 /opt/utils/m05_getInput4matrixQTL.py --pre ${plink_pre}
        Rscript /opt/utils/m06_matrixQTL.R ${plink_pre}.SNP.txt ${plink_pre}.final_cyto.txt ${plink_pre}.final_covariates.txt \
                ${plink_pre}.qtl.txt
        gzip ${plink_pre}.qtl.txt
        """
}


process cell_qtl {
    publishDir "${params.work_path}/f05_cellQTL", mode:"copy"

    input:
        tuple val(pca_pre), path(pca_files)
        tuple val(plink_pre), path(qtl_files)

    output:
        tuple val(plink_pre), path("cell/${plink_pre}*"), emit: cell_qtl_out

    script:
        """
        python3 /opt/utils/m07_cell_qtl.py --pre ${plink_pre}
        """

}


process qq_pplot {
    publishDir "${params.work_path}/f06_QQplot", mode:"copy"

    input:
        tuple val(plink_pre), path(cell_qtl_files)

    output:
        tuple val(plink_pre), path("${plink_pre}*.png")
    
    script:
        """
        python3 /opt/utils/m08_qqplot.py --pre ${plink_pre}
        """
}


process clump {
    publishDir "${params.work_path}/f07_clump", mode:"copy"

    input:
        tuple val(pca_pre), path(pca_files)
        tuple val(plink_pre), path(cell_qtl_files)
        val pval

    output:
        tuple val(plink_pre), path("${plink_pre}*.clumped"), emit: clump_out

    script:
        """
        python3 /opt/utils/m09_clump_qtl.py --pre ${plink_pre} --pval ${pval}
        """
}


process anno_clump {
    publishDir "${params.work_path}/f08_clump_anno", mode:"copy"

    input:
        tuple val(plink_pre), path(cell_qtl_files)
        each path(anno)
        each path(anno_idx)
        val(pval)

    output:
        tuple val(plink_pre), path("${plink_pre}*"), emit: anno_clump_out

    script:
        """
        python3 /opt/utils/m10_annotate_clump_qtl.py --pre ${plink_pre} --pval ${pval}  --anno ${anno}
        """

}



workflow {
    // step1 QC
    Channel
        .fromPath(params.vcf_fn)
        .ifEmpty { exit 1, "Cannot find VCF file : ${params.vcf_fn}"}
        .set {vcf_fn}
        
    Channel
        .fromPath(params.id_map)
        .ifEmpty { exit 1, "Cannot find id map file: ${params.id_map}"}
        .set {id_map_fn}
    
    QC(vcf_fn, id_map_fn)

    // step2 calculate PCA
    Channel
        .fromFilePairs("${params.ld_ref}*.{bed,bim,fam}", size:3)
        .ifEmpty { exit 1, "No matching plink files: ${params.ld_ref}" }
        .set {ld_ref}

    Channel
        .fromPath(params.high_ld)
        .ifEmpty { exit 1, "Cannot find high ld reference file: ${params.high_ld}"}
        .set {high_ld}


    PCA(QC.out.qc_out, ld_ref, high_ld)

    // step 3 get covariation files
    Channel
        .fromPath(params.cyto_fn)
        .ifEmpty { exit 1, "Cannot find cell abundance file: ${params.cyto_fn}"}
        .set {cyto_fn}
    
    get_cova_fn(PCA.out.pca_out, id_map_fn, cyto_fn)

    // step 4 run matrix qtl
    matrix_qtl(QC.out.qc_out, PCA.out.pca_out, get_cova_fn.out.cova_out)

    // step5 split qtl to cell type specific 
    cell_qtl(PCA.out.pca_out, matrix_qtl.out.qtl_out)

    // step 6 generate QQ plot
    qq_pplot(cell_qtl.out.cell_qtl_out)

    // step 7 clump
    clump(PCA.out.pca_out, cell_qtl.out.cell_qtl_out, params.clump_pval)

    // step 8 
    Channel
        .fromPath(params.anno)
        .ifEmpty { exit 1, "Cannot find annotation file: ${params.anno}"}
        .set {anno}
    
    anno_idx_fn = params.anno + '.tbi'
    Channel
        .fromPath(anno_idx_fn)
        .ifEmpty { exit 1, "Cannot find annotation file: ${anno_idx_fn}"}
        .set {anno_idx}

    anno_clump(clump.out.clump_out, anno, anno_idx, params.clump_pval)
}   