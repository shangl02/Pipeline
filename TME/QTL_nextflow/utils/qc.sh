path=/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/
work_path=$path/f02_QC/plink1
plink=$path/software/plink
plink2=$path/software/plink2

#------------------- remove low quality snps and samples -----------
# 1. VCF to plink format 
$plink --vcf $path/tcga.vcf.gz --make-bed --update-sex $work_path/gender.txt --pheno $work_path/pheno.txt --out $work_path/tcga_raw
$plink --bfile $work_path/tcga_raw --missing --out $work_path/f01_missing/tcga_raw_miss
# extract normal samples (10807)
$plink --bfile $work_path/tcga_raw --make-bed --keep $work_path/study_samples.txt --out $work_path/tcga
# 2. genotype call rates and probe call rates
$plink --bfile $work_path/tcga --missing --out $work_path/f01_missing/tcga_miss

# 3.1 bad samples with call rate < 95%
tail -n+2 $work_path/f01_missing/tcga_miss.imiss | awk '{ if ($6>0.05) print $1}' > $work_path/f01_missing/lowCallRate_sample.txt
tail -n+2 $work_path/f01_missing/tcga_miss.imiss | awk '{ if ($6<=0.05) print $1,$1}' > $work_path/f01_missing/highCallRate_sample.txt
# 3.2 SNPs with call rate < 0.95 (49432 snps failed)
tail -n+2 $work_path/f01_missing/tcga_miss.lmiss | awk '{ if ($5>0.05) print $2}' > $work_path/f01_missing/lowCallRate_marker.txt
# 3.3. get duplicate markers (2)
perl /lustre/workspace/projects/ECD/bzhang2/pipeline/GWAS/duplicate_marker.pl $work_path/f01_missing/tcga_miss.lmiss $work_path/tcga.bim > $work_path/f01_missing/duplicate_marker.txt
# 3.4 merge all markers that need to be removed (49434)
cat $work_path/f01_missing/duplicate_marker.txt $work_path/f01_missing/lowCallRate_marker.txt | sort | uniq > $work_path/f01_missing/tobeRemoved_marker.txt
# 3.5 check heterozygosity (10494 samples)
$plink --bfile $work_path/tcga --het --out $work_path/f01_missing/tcga_het
# 3.6 remove call rate (855752 snps, 10494 sample)
$plink --bfile $work_path/tcga --exclude $work_path/f01_missing/tobeRemoved_marker.txt --keep $work_path/f01_missing/highCallRateHet_sample.txt --make-bed --out $work_path/f01_missing/tcga_rmSp_rmSNP

# 4. check sex (855,752, 10458 samples pass the test)
$plink --bfile $work_path/f01_missing/tcga_rmSp_rmSNP --check-sex 0.40 0.4 --out $work_path/f02_checksex/tcga_rmSp_rmSNP_sexcheck

# tail -n+2 $work_path/f02_checksex/tcga_rmSp_rmSNP_sexcheck.sexcheck | awk '{ if ($5 == "OK") print $1,$1}' > $work_path/f02_checksex/keep_sex.txt

$plink --bfile $work_path/f01_missing/tcga_rmSp_rmSNP --keep $work_path/f02_checksex/keep_sex.txt --make-bed --out $work_path/f02_checksex/tcga_rmSp_rmSNP_sex 


#####----------- merge with 1000 genomics and do pca----------------

# Extract autosomal snps (820,890 snps, 10458 samples)
$plink --bfile $work_path/f02_checksex/tcga_rmSp_rmSNP_sex --chr 1-22 --make-bed --out $work_path/f03_1kg/plinkfiletorun.autosomal

fgrep rs $work_path/f03_1kg/plinkfiletorun.autosomal.bim | awk '{print $2}' - > $work_path/f03_1kg/plinkfiletorun.rsid_names.txt
# Extract 1kg autosomal snps (793,686 snps, 2504 samples)
$plink --bfile /hpc/grid/hgcb/workspace/projects/P002_reference_information/plink/1000G_6PoP_MAF0.05 --extract $work_path/f03_1kg/plinkfiletorun.rsid_names.txt --make-bed --out $work_path/f03_1kg/1kg_phase1_all.rsids.autosomal

# Obtain SNPs present in both files (793,686 snps)
awk '{print $2}' $work_path/f03_1kg/1kg_phase1_all.rsids.autosomal.bim > $work_path/f03_1kg/1kg_phase1_all.rsids_names.txt

# Extract 1KG SNPs from root (793, 686 snps, 10458)
$plink --bfile $work_path/f03_1kg/plinkfiletorun.autosomal --extract $work_path/f03_1kg/1kg_phase1_all.rsids_names.txt --make-bed --out $work_path/f03_1kg/plinkfiletorun.intersection

# Dry run bmerge to identify SNPs PLINK will fail on (12 variants with multi allele)
$plink --bfile $work_path/f03_1kg/plinkfiletorun.intersection --bmerge $work_path/f03_1kg/1kg_phase1_all.rsids.autosomal --merge-mode 6 --out $work_path/f03_1kg/plinkfiletorun.1KG.failures

# Add variants with multiple positions to missnp (4 snps wiht multi position)
fgrep \'rs $work_path/f03_1kg/plinkfiletorun.1KG.failures.log |\
awk '{print $7}' |\
sed -e "s/'//g" -e "s/\.//g" > $work_path/f03_1kg/plinkfiletorun.1KG.failures.multiple.positions.txt
# merge multi-allel and multi-position (16)
cat $work_path/f03_1kg/plinkfiletorun.1KG.failures.missnp $work_path/f03_1kg/plinkfiletorun.1KG.failures.multiple.positions.txt > $work_path/f03_1kg/plinkfiletorun.1KG.failures.multiple.positions.missnp  

# Exclude mismatched SNPs and variants with multiple positions and alleles (793, 671 snps, 10458 samples)
$plink --bfile $work_path/f03_1kg/plinkfiletorun.intersection --exclude $work_path/f03_1kg/plinkfiletorun.1KG.failures.multiple.positions.missnp --make-bed --out $work_path/f03_1kg/plinkfiletorun.intersection_for_merge


# Merge root and 1KG (793,686 snps, 12962 samples)  
$plink --bfile $work_path/f03_1kg/plinkfiletorun.intersection_for_merge --bmerge $work_path/f03_1kg/1kg_phase1_all.rsids.autosomal --out $work_path/f03_1kg/plinkfiletorun.1kg.pop_strat

# Filter missing variants, rare variants and HWE , exclude high LD region (312,916 snps, 12962 samples)
$plink --bfile $work_path/f03_1kg/plinkfiletorun.1kg.pop_strat --geno 0.01 --maf 0.01 --hwe 0.000001 --exclude range /lustre/workspace/projects/ECD/bzhang2/data/high-LD-regions.txt --make-bed --out $work_path/f03_1kg/plinkfiletorun.1kg.pop_strat.for_prune

# LD Pruning (66,688 SNPs, 12962 samples)
$plink2 --bfile $work_path/f03_1kg/plinkfiletorun.1kg.pop_strat.for_prune --indep-pairwise 1500 150 0.2 --out $work_path/f03_1kg/plinkfiletorun.1kg.pop_strat.prune

$plink2 --bfile $work_path/f03_1kg/plinkfiletorun.1kg.pop_strat.for_prune --extract $work_path/f03_1kg/plinkfiletorun.1kg.pop_strat.prune.prune.in --make-bed --out $work_path/f03_1kg/plinkfiletorun.1kg.LD_pop_strat


### extract 1kg for flashpca
mkdir $work_path/f04_pca

$plink2 --bfile $work_path/f03_1kg/plinkfiletorun.1kg.LD_pop_strat --chr 1-22 --keep /hpc/grid/hgcb/workspace/projects/P002_reference_information/plink/1000G_6PoP_MAF0.05.fam --make-bed --out $work_path/f04_pca/1kg_Precision

### extract study population (10458 samples,66688 snps)
$plink2 --bfile $work_path/f03_1kg/plinkfiletorun.1kg.LD_pop_strat --chr 1-22 --keep $work_path/f03_1kg/plinkfiletorun.intersection.fam --make-bed --out $work_path/f04_pca/tcga_Pruned

/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/software/flashpca_x86-64 --bfile $work_path/f04_pca/1kg_Precision --suffix _1kg --outload $work_path/f04_pca/1kg_loadings.txt --outmeansd $work_path/f04_pca/1kg_meansd.txt
### project to study population 
/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/software/flashpca_x86-64 --bfile $work_path/f04_pca/tcga_Pruned --project --inmeansd $work_path/f04_pca/1kg_meansd.txt --outproj $work_path/f04_pca/tcga_projections.txt --inload $work_path/f04_pca/1kg_loadings.txt --suffix _study -v

### extract white race samples


#### relatedness using non-pruned QCed plink file 
bsub -q medium -o log.king.relate -n 16 -R "span[ptile=16]" "/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/software/king -b $work_path/f01_missing/tcga_rmSp_rmSNP.bed --prefix tcga.qced --related --degree 4 --cpus 16 --rplot"

bsub -q medium -o log.king.cluster -n 16 -R "span[ptile=16]" "/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/software/king -b $work_path/f01_missing/tcga_rmSp_rmSNP.bed --prefix tcga.qced --cluster --degree 4 --cpus 16"

bsub -q medium -o log.king.build -n 16 -R "span[ptile=16]" "/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/software/king -b $work_path/f01_missing/tcga_rmSp_rmSNP.bed --prefix tcga.qced --build --degree 3 --cpus 16"



#=============   prepare imputation ===============================
$plink --freq --bfile $work_path/f02_checksex/tcga_rmSp_rmSNP_sex --out $work_path/pre_impute_check/tcga_impute_freq

perl $work_path/pre_impute_check/HRC-1000G-check-bim.pl -b $work_path/f02_checksex/tcga_rmSp_rmSNP_sex.bim -f $work_path/pre_impute_check/tcga_impute_freq.frq -r $work_path/pre_impute_check/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h sh Run-plink.sh

# Quality Control for HRC, 1000G and CAAPA imputation
$plink --bfile $work_path/f02_checksex/tcga_rmSp_rmSNP_sex --exclude $work_path/pre_impute_check/Exclude-tcga_rmSp_rmSNP_sex-HRC.txt --make-bed --out $work_path/pre_impute_check/TEMP1
$plink --bfile $work_path/pre_impute_check/TEMP1 --update-map $work_path/pre_impute_check/Chromosome-tcga_rmSp_rmSNP_sex-HRC.txt --update-chr --make-bed --out $work_path/pre_impute_check/TEMP2
$plink --bfile $work_path/pre_impute_check/TEMP2 --update-map $work_path/pre_impute_check/Position-tcga_rmSp_rmSNP_sex-HRC.txt --make-bed --out $work_path/pre_impute_check/TEMP3
$plink --bfile $work_path/pre_impute_check/TEMP3 --flip $work_path/pre_impute_check/Strand-Flip-tcga_rmSp_rmSNP_sex-HRC.txt --make-bed --out $work_path/pre_impute_check/TEMP4
$plink --bfile $work_path/pre_impute_check/TEMP4 --reference-allele $work_path/pre_impute_check/Force-Allele1-tcga_rmSp_rmSNP_sex-HRC.txt --make-bed --out $work_path/pre_impute_check/tcga_rmSp_rmSNP_sex-updated

# transfer plink to vcf, which will be submitted to Michigan Imputation
for i in {1..22}; do $plink --bfile $work_path/pre_impute_check/tcga_rmSp_rmSNP_sex-updated --a2-allele $work_path/pre_impute_check/Force-Allele1-tcga_rmSp_rmSNP_sex-HRC.txt --recode vcf-iid bgz --chr $i --out $work_path/pre_impute_check/chr/tcga_rmSp_rmSNP_sex_updated_chr$i;done;

$plink --bfile $work_path/pre_impute_check/tcga_rmSp_rmSNP_sex-updated --a2-allele $work_path/pre_impute_check/Force-Allele1-tcga_rmSp_rmSNP_sex-HRC.txt --keep $work_path/f02_checksex/males.txt --recode vcf-iid bgz --chr 23 --out $work_path/pre_impute_check/chr/tcga_rmSp_rmSNP_sex_updated_chr23_male
$plink --bfile $work_path/pre_impute_check/tcga_rmSp_rmSNP_sex-updated --a2-allele $work_path/pre_impute_check/Force-Allele1-tcga_rmSp_rmSNP_sex-HRC.txt --keep $work_path/f02_checksex/females.txt --recode vcf-iid bgz --chr 23 --out $work_path/pre_impute_check/chr/tcga_rmSp_rmSNP_sex_updated_chr23_female
gunzip -c tcga_rmSp_rmSNP_sex_updated_chr23_male.vcf.gz | sed 's/^23/X/g' - | bgzip  > tcga_rmSp_rmSNP_sex_updated_chrX_male.vcf.gz
gunzip -c tcga_rmSp_rmSNP_sex_updated_chr23_female.vcf.gz | sed 's/^23/X/g' - | bgzip  > tcga_rmSp_rmSNP_sex_updated_chrX_female.vcf.gz
# check vcf to make sure the allele type match the reference genome
for i in {22..23}; do python2 checkVCF.py -r /lustre/scratch/lis262/hs37d5.fa -o chr$i /hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f02_QC/plink1/pre_impute_check/chr/tcga_rmSp_rmSNP_sex_updated_chr$i.vcf.gz;done;

python2 checkVCF.py -r /lustre/scratch/lis262/hs37d5.fa -o test /hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f02_QC/plink1/pre_impute_check/chr/tcga_rmSp_rmSNP_sex_updated_chr22.vcf.gz


#================ the following is for preparing files for matrixQTL
$plink --bfile /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink/merge/autosomal --extract /hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f03_gwas/test/id.txt --make-bed --out /hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f03_gwas/test/test

$plink --bfile /hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f02_QC/plink1/pre_impute_check/tcga_rmSp_rmSNP_sex-updated --extract /hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f03_gwas/test/id.txt --make-bed --out /hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f03_gwas/test/test

#======= change vcf to plink format
# extract interger genotype
for i in {21..22}; do $plink2 --vcf /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/chr$i.dose.vcf.gz --extract-if-info "R2 >= 0.3" --recode A-transpose --out /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink/chr$i --recode-allele /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/chr$i.ale.tsv && bgzip /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink/chr$i.traw;done;
# extract dosage genotype
for i in {1..22}; do $plink2 --vcf /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/chr$i.dose.vcf.gz --extract-if-info "R2 >= 0.3" --recode A-transpose --out /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink/chr$i --recode-allele /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/chr$i.ale.tsv && gzip /hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink/chr$i.traw;done;

#=========== vcf to plink2 format for case control analysis

$plink2 --pfile plink/chr22.dose --glm --const-fid 0 --pheno covari.txt --pheno-name pheno  --covar-name PC1,PC2,PC3,PC4,gender,age --covar-variance-standardize --threads 4 --out chr22