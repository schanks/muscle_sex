library(data.table)


inverse.normalize <- function(i) {
  stopifnot(is.vector(i, mode = "numeric"))
  qnorm((rank(i,
    na.last = "keep",
    ties.method = "random"
  ) - 0.5) / sum(!is.na(i)))
}



results=fread("../all.expr.results")
means=fread("../cell_means_bygene.tab")
map=fread("/net/dumbo/home/dciotlos/R_Scripts/ensembl_to_entrez_V1_apr21.csv")

results$ENSEMBL=unlist(strsplit(results$gene, "[.]",))[seq(1, nrow(results)*2, by=2)]
results=merge(results, map[,c("ENSEMBL","Entrez_ID")], by="ENSEMBL")
results=merge(results, means, by=c("gene","cell"))

#All genes - invnorm -log10 p-values
all=results[,c("Entrez_ID","log2FoldChange","pvalue", "meanCount","cell")]

t1_all=all[which(all$cell=="Type_1"), c("Entrez_ID","pvalue","log2FoldChange","meanCount")]
t2a_all=all[which(all$cell=="Type_2a"), c("Entrez_ID","pvalue","log2FoldChange","meanCount")]
t2x_all=all[which(all$cell=="Type_2x"), c("Entrez_ID","pvalue","log2FoldChange","meanCount")]

write.table(t1_all[,c("Entrez_ID","pvalue","log2FoldChange","meanCount")], "~/t1_input_pval.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t2x_all[,c("Entrez_ID","pvalue","log2FoldChange","meanCount")], "~/t2x_input_pval.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t2a_all[,c("Entrez_ID","pvalue","log2FoldChange","meanCount")], "~/t2a_input_pval.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

