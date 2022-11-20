# Created by Shangzhong.Li@pfizer.com on 2021/02/25
# run regression for peer
args = commandArgs(trailingOnly=TRUE)


tpm_fn = args[1]  # rows are genes, columns are samples
cyto_fn = args[2] # rows are samples, columns are cells
out_fn = args[3]

tpm_df = read.table(tpm_fn,sep='\t',header=T,row.names=1)
ind = (rowSums(tpm_df>0.1)>=10)
inv.normal = function(x) qnorm((rank(x)-0.5)/length(x))
invnorm = apply(tpm_df[ind, ], 1, inv.normal)
invnorm <- invnorm[sort(row.names(invnorm)),]
genes = colnames(invnorm)
genes = gsub("-", "\\.", genes)
colnames(invnorm) = genes
invnorm = as.data.frame(invnorm)
duplicated <- duplicated(row.names(invnorm))
invnorm = invnorm[!duplicated,]
#----- cytoreaons ----
# cyto_fn = '/home/rstudio/p001_RNASeq/tme/cyto.txt'
cyto_df = read.table(cyto_fn,sep='\t',header=T,row.names=1)
cells = colnames(cyto_df)
#===== merge data and run glm =========
merge_df = merge(invnorm, cyto_df, by="row.names")
result = c()
for (g in genes) {
  ivnames <- paste(cells,sep="",collapse="+")
  form <- formula(paste(g,"~",ivnames))
  tryCatch(
    {md <- glm(form, data=merge_df)
    residue <- md$residuals
    result <- rbind(result, residue)},
    error=function(e){cat(g,"ERROR:",conditionMessage(e),"\n")}
  )
}
rownames(result) = genes
colnames(result) = merge_df$Row.names
# out_fn = '/home/rstudio/p001_RNASeq/tme/residue.tsv'
write.table(result,out_fn,sep="\t",quote=F)
