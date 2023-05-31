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



include {saige_step1_merge_gt; saige_step2_bgen} from './utils_saige'

workflow {
    Channel
        .fromFilePairs("${params.step1_gt_pre}*.{bed,bim,fam}", size:3)
        .ifEmpty { error "No matching plink files" }
        .set {pfiles}

    Channel.fromPath(params.pheno_cova_fn)
        .ifEmpty { exit 1, "Cannot find pheno_cova file : ${params.pheno_cova_fn}" } 
        .set {pheno_fn}

    saige_step1_merge_gt(pfiles, pheno_fn, params.covariates, params.cate_cova, params.trait_type, ava_cpu)


    Channel
        .fromFilePairs("${params.step2_gt_pre}*.{bgen,bgen.bgi,sample}", flat : true, size:3)
        .ifEmpty {error "No mathcing bgen files"}
        .set {bgenFiles}


    saige_step2_bgen(saige_step1_merge_gt.out.step1_out, bgenFiles)

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Pipeline failed" )
}