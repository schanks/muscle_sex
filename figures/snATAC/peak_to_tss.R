library(data.table)

celltypes=c("Type_1","Type_2a","Type_2x","Endothelial","Mesenchymal_Stem_Cell","Macrophage","Neuromuscular_junction","Neuronal","Satellite_Cell","Smooth_Muscle")

window=1000

for (c1 in celltypes){
	genes=fread(paste("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/final_drop10nuc/results/",c1,".SEX.M.results.tab", sep=""))
	genes$N_peak=rep(0, nrow(genes))
	genes$F_peak=rep(0, nrow(genes))
	genes$M_peak=rep(0, nrow(genes))
	genes$N_sig_peak=rep(0, nrow(genes))
	genes$F_sig_peak=rep(0, nrow(genes))
	genes$M_sig_peak=rep(0, nrow(genes))
	info=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/gencode.v30.gene_lengths_tss.tsv")
	colnames(info)[4]="gene"
	print(dim(genes))
	genes=merge(genes, info[,c("gene","tss_start")], by="gene")
	print(dim(genes))
	atac=fread(paste("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.ATAC/final_drop10nuc/results/",c1,".SEX.M.results.tab", sep=""))
	atac$fdr=p.adjust(atac$pvalue, method="fdr")
	atac$chrom=unlist(strsplit(atac$peak, ":"))[seq(1, nrow(atac)*3, by=3)]
	atac$start=as.numeric(unlist(strsplit(atac$peak, ":"))[seq(2, nrow(atac)*3, by=3)])
	atac$stop=as.numeric(unlist(strsplit(atac$peak, ":"))[seq(3, nrow(atac)*3, by=3)])
	atac$midpt=(atac$stop+atac$start)/2
	atac$direction=as.numeric(atac$log2FoldChange>0)
	atac$sig=as.numeric(atac$fdr<0.05)

	for (i in 1:nrow(genes)){
		if (genes$gene_strand[i]=="+"){
			peaks=atac[which(atac$chrom==genes$chrom[i] & genes$tss_start[i]-atac$midpt>0 & genes$tss_start[i]-atac$midpt<=window),]
		}
		else{
			peaks=atac[which(atac$chrom==genes$chrom[i] & atac$midpt-genes$tss_start[i]>0 & atac$midpt-genes$tss_start[i]<=window),]
		}
	genes$N_peak[i]=nrow(peaks)
	genes$F_peak[i]=sum(peaks$direction==0)
	genes$M_peak[i]=sum(peaks$direction)
	genes$N_sig_peak[i]=sum(peaks$sig)
	genes$F_sig_peak[i]=sum(peaks$sig==1 & peaks$direction==0)
	genes$M_sig_peak[i]=sum(peaks$sig==1 & peaks$direction==1)
	Sys.sleep(0.01)
	print(i/nrow(genes))
	}

	write.table(genes, paste(c1, ".1kb.peak.results", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}












