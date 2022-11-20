
# Created by Shangzhong.Li@pfizer.com on 2021/02/25

args = commandArgs(trailingOnly=TRUE)

library(peer)

residue_fn = args[1]
residue_df = read.table(residue_fn,sep="\t",header=T,row.names=1,quote="")

out_pre = args[2]

WriteTable <- function(data, filename, index.name) {
  datafile <- file(filename, open = "wt")
  on.exit(close(datafile))
  header <- c(index.name, colnames(data))
  writeLines(paste0(header, collapse="\t"), con=datafile, sep="\n")
  write.table(data, datafile, sep="\t", col.names=FALSE, quote=FALSE)
}


inv.normal <- function(x) qnorm((rank(x)-0.5)/length(x))
invnorm <- apply(residue_df, 1, inv.normal)

N <- nrow(invnorm)

k <- 15 * (1 + min(1, floor(N/150)) + min(1, floor(N/250)) + min(1, floor(N/350)))
model <- PEER()
invisible(PEER_setNk(model, k))
invisible(PEER_setPhenoMean(model, invnorm))
invisible(PEER_setPriorAlpha(model, 0.001, 0.01))
invisible(PEER_setPriorEps(model, 0.1, 10.0))
invisible(PEER_setNmax_iterations(model, 1000))

PEER_update(model)

X <- PEER_getX(model)  # samples x PEER factors
A <- PEER_getAlpha(model)  # PEER factors x 1
R <- t(PEER_getResiduals(model))  # genes x sample

cols <- paste0("PEER",1:ncol(X))
rownames(X) <- rownames(invnorm)
colnames(X) <- cols
rownames(A) <- cols
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha
rownames(R) <- colnames(invnorm)
colnames(R) <- rownames(invnorm)



WriteTable(t(X), paste0(out_pre, ".PEER_covariates.txt"), "ID")  # format(X, digits=6)
WriteTable(A, paste0(out_pre, ".PEER_alpha.txt"), "ID")
WriteTable(R, paste0(out_pre, ".PEER_residuals.txt"), "ID")