# vcf_fn=/lustre/workspace/projects/BLCA/Results/WES/vcf/anno/wes_germ.norm.rsid.hg38.vcf.gz
vcf_fn=/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/imp_WES/germ_newSp38_imp_rsid.vcf.gz
work_path=/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/imp_WES
code_path=~/Code/Pipeline/TME/QTL_pipeline/utils/
id_map=$work_path/wes_id_map.tsv
cyto_fn=/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/WGS/Javline_cytoreason_table.tsv
anno=/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/imp_WES/anno/imp_WES.vep.tsv.gz
# step1
# python $code_path/m01_process_gt.py -i $vcf_fn -p $work_path
# step2
# python $code_path/m02_clean_samples.py -p $work_path --id_map $id_map
# step3
# python $code_path/m03_calculate_PC.py -p $work_path
# step4
# python $code_path/m04_get_covariate_file.py -p $work_path --id_map $id_map -d $cyto_fn
# step5
# python $code_path/m05_getInput4matrixQTL.py -p $work_path
# step6
# python $code_path/m06_run_matrixQTL.py -p $work_path
# step7
# python $code_path/m07_cell_qtl.py -p $work_path
# step8 
# python $code_path/m08_qqplot.py -p $work_path
# step9
python $code_path/m09_clump_qtl.py -p $work_path --pval 5e-7
# step10
python $code_path/m10_annotate_clump_qtl.py -p $work_path -a $anno


