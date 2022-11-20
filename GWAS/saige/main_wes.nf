#!/usr/bin/env nextflow

/* Enable DSL2 syntax */
nextflow.enable.dsl=2
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


include {get_chrom_bfiles; get_chrom_vcf_files; saige_step1_chrom_gt; saige_step2_vcf} from './utils_saige'



workflow {    
    Channel
        .of(1..22)
        .set {chrom}
    
    // Step1
    chrom
        .map {it -> [it, "${params.step1_gt_pre}${it}"]}
        .set {chrom_pfile_pre}
    
    get_chrom_bfiles(chrom_pfile_pre)
        .splitCsv(by:5, sep:'\n')
        .map {it.flatten().toList()}
        .set {chrom_bfiles}
        
    Channel.fromPath(params.pheno_cova_fn)
        .ifEmpty { exit 1, "Cannot find pheno_cova file : ${params.pheno_cova_fn}" } 
        .set {pheno_fn}

    saige_step1_chrom_gt(chrom_bfiles, pheno_fn, params.covariates, params.cate_cova, params.trait_type, ava_cpu)
    
    // Step 2
    chrom
        .map {it -> [it, "${params.step2_gt_pre}${it}"]}
        .set {chrom_vcf_pre}
    
    get_chrom_vcf_files(chrom_vcf_pre)
        .splitCsv(by:4, sep:'\n')
        .map {it.flatten().toList()}
        .set {chrom_vcf_files}
    
    saige_step1_chrom_gt.out.step1_out
        .join(chrom_vcf_files, by:0)
        .set {input4step2}

    saige_step2_vcf(input4step2, params.vcf_field)

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Pipeline failed" )
}