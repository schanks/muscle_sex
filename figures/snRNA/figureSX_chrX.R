library(data.table)
library(ggplot2)

sig_count<-function(x){
	return(sum(x<0.05))
}

celltypes=c("Type_1","Type_2a","Type_2x","Mesenchymal_Stem_Cell","Satellite_Cell","Neuromuscular_junction","Endothelial","Smooth_Muscle","Macrophage","Neuronal")

cell_map=as.data.frame(celltypes)
colnames(cell_map)="cell"
cell_map$list=c(10,11,12,4,7,5,2,8,3,6)
cell_map$mean=rep(0, nrow(cell_map))
cell_map$total=rep(0, nrow(cell_map))

eds.orig<-readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/data/eds.orig.Rds")
pds.orig<-readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/data/pds.orig.Rds")

for (i in 1:nrow(cell_map)){
	index=cell_map$list[i]
	pds=pds.orig[[index]]
	pds=pds[which(pds$tot_nuc>10),]
	samps.drop10=pds$labelcode
	
	eds=eds.orig[[index]]
	gene=as.data.frame(eds$gene)
	colnames(gene)=c("gene")
	count.tmp=eds[,-1]
	nsamp=ncol(count.tmp)
	ct.non0=rowSums(count.tmp!=0)
	gene$prop.non0=ct.non0/nsamp
	gene=gene[which(gene$prop.non0>=0.25),]
	eds=eds[which(is.element(eds$gene, gene$gene)),]

	eds=eds[,which(is.element(colnames(eds), samps.drop10))]
	cell_map$mean[i]=mean(colSums(eds))
	cell_map$total[i]=sum(colSums(eds))
}

results=fread("/net/snowwhite/home/schanks/muscle/snrna/all.expr.results")
results=results[which(results$fdr<0.05),]
results=results[which(results$chrom!="chrY"),]
results$CHR=rep("Autosomal", nrow(results))
results$CHR[which(results$chrom=="chrX")]="X"

cell_tot=aggregate(results$fdr, by=list(results$cell), FUN=length)
colnames(cell_tot)=c("cell","all")
cell_tot$female=aggregate(results$fdr[which(results$log2FoldChange<0)], by=list(results$cell[which(results$log2FoldChange<0)]), FUN=length)$x
cell_tot$male=aggregate(results$fdr[which(results$log2FoldChange>0)], by=list(results$cell[which(results$log2FoldChange>0)]), FUN=length)$x

cell_tot$faut=aggregate(results$fdr[which(results$log2FoldChange<0 & results)])
cell_tot$maut

cell_tot$fx
cell_tot$mx

cell_tot$propxf=cell_tot$fx/cell_tot$female
cell_tot$propxm=cell_tot$mx/cell_tot$male


cell_tot$sigtot=aggregate(results$fdr, by=list(results$cell), FUN=sig_count)$x
cell_tot$proptot=cell_tot$sigtot/cell_tot$tested

cell_tot=merge(cell_tot, cell_map[,c("mean","total","cell")], by=c("cell"))









library(data.table)
library(ggplot2)
library(gridExtra)

sig_count<-function(x){
	return(sum(x<0.05))
}

results=fread("/net/snowwhite/home/schanks/muscle/snrna/all.expr.results")
results=results[which(results$fdr<0.05),]
results=results[which(results$chrom!="chrY")]
results$CHR=rep("Autosomal", nrow(results))
results$CHR[which(results$chrom=="chrX")]="X"

results$Direction=rep("Female-biased", nrow(results))
results$Direction[which(results$log2FoldChange>0)]="Male-biased"
results$cell[which(results$cell=="Mesenchymal_Stem_Cell")]="FAP"
results$cell[which(results$cell=="Neuromuscular_junction")]="NMJ"
results$cell=gsub("_"," ", results$cell)

tiff("~/plot.tiff", height=180, width=180, units="mm", res=300)
ggplot(results, aes(x=Direction))+geom_bar(position="fill", aes(fill=CHR))+facet_wrap(cell~.)+scale_fill_manual(values=c("#e41a1c","#377eb8","gray"))+theme_bw()+theme(axis.title.x=element_blank(), strip.background=element_rect(fill=NA, color=NA))+ylab("Proportion of tested genes\nby differential expression status")
dev.off()



