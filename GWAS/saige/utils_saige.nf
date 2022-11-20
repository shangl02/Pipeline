
process saige_step1_merge_gt {
    tag "${pfile_prefix}"
    publishDir "${params.outdir}/step1", mode: "copy"

    input:
        tuple val(pfile_prefix), path(pfiles)
        each path(pheno_cova_fn)
        val covar
        val cate_cova
        val trait_type
        val thread
    
    output:
        tuple path("${pfile_prefix}_step1.rda"), path("${pfile_prefix}_step1.varianceRatio.txt"), emit: step1_out

    script:
        if(params.trait_type == "quantitative")
        """
        step1_fitNULLGLMM.R \
            --plinkFile=${pfile_prefix} \
            --phenoFile=${pheno_cova_fn} \
            --phenoCol=pheno \
            --covarColList=${covar} \
            --qCovarColList=${cate_cova} \
            --sampleIDColinphenoFile=IID \
            --traitType=${trait_type} \
            --invNormalize=TRUE \
            --outputPrefix=${pfile_prefix}_step1 \
            --nThreads=${thread} \
            --memoryChunk=2 \
            --LOCO=FALSE  \
            --tauInit=1,0 \
            --IsOverwriteVarianceRatioFile=TRUE
        """
    else if(params.trait_type=="binary")
        """
        step1_fitNULLGLMM.R \
            --plinkFile=${pfile_prefix} \
            --phenoFile=${pheno_cova_fn} \
            --phenoCol=pheno \
            --covarColList=${covar} \
            --qCovarColList=${cate_cova} \
            --sampleIDColinphenoFile=IID \
            --traitType=${trait_type} \
            --invNormalize=FALSE \
            --outputPrefix=${pfile_prefix}_step1 \
            --nThreads=9 \
            --memoryChunk=2 \
            --LOCO=FALSE  \
            --tauInit=0,0 \
            --IsOverwriteVarianceRatioFile=TRUE
        """
}


process saige_step2_bgen {
    tag "${bgen_pre}"
    publishDir "${params.outdir}/step2", mode: 'copy'

    input:
        each path(rda_vari)
        tuple val(bgen_pre), path(bgen_file), path(bgen_idx), path(bgen_sample)

    output:
        val(bgen_pre), emit: step2_out

    script:
        """
        step2_SPAtests.R \
            --bgenFile=${bgen_file}    \
            --bgenFileIndex=${bgen_idx}   \
            --sampleFile=${bgen_sample} \
            --minMAF=0 \
            --minMAC=3 \
            --GMMATmodelFile=${rda_vari[0]} \
            --varianceRatioFile=${rda_vari[1]}   \
            --SAIGEOutputFile=${bgen_pre}_asso \
            --LOCO=FALSE  # --numLinesOutput=1000
        """
}


