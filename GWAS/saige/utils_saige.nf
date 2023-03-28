
// get chromosome wise plink files
process get_chrom_bfiles {
    input:
        tuple val(chrom), val(prefix)
    
    output:
        stdout

    script:
        """
        echo ${chrom}
        basename ${prefix}
        echo ${prefix}.bed
        echo ${prefix}.bim
        echo ${prefix}.fam
        """
}

// get chromosome wise vcf files
process get_chrom_vcf_files {
    input:
        tuple val(chrom), val(prefix)
    
    output:
        stdout

    script:
        """
        echo ${chrom}
        basename ${prefix}
        echo ${prefix}.vcf.gz
        echo ${prefix}.vcf.gz.csi
        """
}

// saige step1, run for each chromosome, use chromosome as index of the input files
process saige_step1_chrom_gt {
    tag "${prefix}"
    publishDir "${params.outdir}/step1", mode: "copy"

    input:
        tuple val(chrom), val(prefix), path(bed), path(bim), path(fam)
        each path(pheno_cova_fn)
        val covar
        val cate_cova
        val trait_type
        val thread

    output:
        tuple val(chrom), path("${prefix}_step1.rda"), path("${prefix}_step1.varianceRatio.txt"), emit: step1_out
    
    script:
        if(params.trait_type == "quantitative")
        """
        step1_fitNULLGLMM.R \
            --plinkFile=${prefix} \
            --phenoFile=${pheno_cova_fn} \
            --phenoCol=pheno \
            --covarColList=${covar} \
            --qCovarColList=${cate_cova} \
            --sampleIDColinphenoFile=IID \
            --traitType=${trait_type} \
            --invNormalize=TRUE \
            --outputPrefix=${prefix}_step1 \
            --nThreads=${thread} \
            --memoryChunk=2 \
            --LOCO=FALSE  \
            --tauInit=1,0 \
            --IsOverwriteVarianceRatioFile=TRUE
        """
        else if(params.trait_type == "binary")
        """
        step1_fitNULLGLMM.R \
            --plinkFile=${prefix} \
            --phenoFile=${pheno_cova_fn} \
            --phenoCol=pheno \
            --covarColList=${covar} \
            --qCovarColList=${cate_cova} \
            --sampleIDColinphenoFile=IID \
            --traitType=${trait_type} \
            --invNormalize=FALSE \
            --outputPrefix=${prefix}_step1 \
            --nThreads=${thread} \
            --memoryChunk=2 \
            --LOCO=FALSE  \
            --tauInit=0,0 \
            --IsOverwriteVarianceRatioFile=TRUE
        """
}


// saige step1, run for a single plink file which merges all variants together.
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

// saige step2 for vcf files, use chromosome as input index
process saige_step2_vcf {
    tag "${vcf_pre}"
    publishDir "${params.outdir}/step2", mode: 'copy'

    input:
        tuple val(chrom), path(rda), path(variation_ratio), val(vcf_pre), path(vcf_file), path(vcf_idx)
        val vcf_field

    output:
        tuple val(vcf_pre), path("${vcf_pre}_asso.txt"), path("${vcf_pre}_asso.txt.index"), emit: step2_out

    script:
        """
        step2_SPAtests.R \
            --vcfFile=${vcf_file}    \
            --vcfFileIndex=${vcf_idx}   \
            --vcfField=${vcf_field} \
            --chrom=${chrom} \
            --minMAF=0 \
            --minMAC=0.5 \
            --GMMATmodelFile=${rda} \
            --varianceRatioFile=${variation_ratio}      \
            --SAIGEOutputFile=${vcf_pre}_asso.txt \
            --LOCO=FALSE               
        """
}

// saige step2 for bgen files, automatically capture all files, no need to use chromosome as index
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


