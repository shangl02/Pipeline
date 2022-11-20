# Created by Shangzhong.Li@pfizer.com on 2021/07/21
# This file run susieR


args = commandArgs(trailingOnly=TRUE)
gwas_fn = args[1]   # gwas summary stats file
ld_fn = args[2]     # LD matrix file for the analysing variants
plot = args[3]      # plot shows the finemap results
out_fn = args[4]    # file that stores the susie results
stat = args[5]     # column name in gwas results indicating the Z score of the variant

library(susieR)
gwas_df = read.csv(gwas_fn,sep='\t')
ld_df = read.table(ld_fn,sep='')

idx_na = which(is.na(gwas_df[,c(stat)]))
if (!is.null(idx_na) & length(idx_na) > 0) {
	gwas_df = gwas_df[-idx_na,]
	ld_df = ld_df[-idx_na, -idx_na]
	print(dim(ld_df))
}

idx_to_rm = c()
while (sum(is.na(ld_df)) > 0) {
  na_count = apply(ld_df,MARGIN=1,function(x) sum(is.na(x)))
  idx_to_remove = which.max(na_count)
  ld_df = ld_df[-idx_to_remove,-idx_to_remove]
  gwas_df = gwas_df[-idx_to_remove,]  
}

R = as.matrix(ld_df)
gz1 = gzfile(gwas_fn,"w")
write.table(gwas_df,gz1,sep='\t',quote=F,row.names=F)
close(gz1)
gz1 = gzfile(ld_fn,"w")
write.table(R,gz1,sep='\t',quote=F,row.names=F,col.names=F)
close(gz1)



fitted_rss = susie_rss(gwas_df[,c(stat)], R, L=10)
cs = summary(fitted_rss)$cs
write.table(cs,out_fn,sep='\t',quote=F,row.names=F)

jpeg(plot)
susie_plot(fitted_rss, y='PIP')
dev.off()