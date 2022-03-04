args <- commandArgs(trailingOnly = TRUE)

path <- args[1]
result <- args[2]

library("MatrixEQTL")
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR
setwd(path)
# setwd('/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/f03_gwas/STAD')
dir = getwd()
SNP_file_name = paste(dir, "/SNP.txt", sep="")
expression_file_name = paste(dir, "/cyto.txt", sep="")
covariates_file_name = paste(dir, "/covariates.txt", sep="")
output_file_name = paste(dir, "/",result, sep="")
pvOutputThreshold = 1
errorCovariance = numeric()
# read snps
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 60,000 rows
snps$LoadFile( SNP_file_name );

maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)
cat('SNPs before filtering:',nrow(snps))
snps$RowReorder(maf>0.001)
cat('SNPs before filtering:',nrow(snps))
# read expression
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$LoadFile(expression_file_name)
# read covariate
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile(covariates_file_name)


me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

# qtl_df = read.csv(output_file_name,sep='\t')
# qqnorm(-10*log10(qtl_df$p.value), pch=1,frame=FALSE)
# qqplot(qtl_df$p.value)
# source('/hpc/grid/hgcb/workspace/projects/GeneticsTools_BinariesAndScripts/R_scripts/qqplotpvaluesII.R')
# qqplot.pvalues(qtl_df$p.value, top.n=500)

