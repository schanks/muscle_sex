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

cell_tot=aggregate(results$fdr, by=list(results$cell), FUN=length)
colnames(cell_tot)=c("cell","tested")
cell_tot$sigtot=aggregate(results$fdr, by=list(results$cell), FUN=sig_count)$x
cell_tot$proptot=cell_tot$sigtot/cell_tot$tested

cell_tot=merge(cell_tot, cell_map[,c("mean","total","cell")], by=c("cell"))

cor(cell_tot$total, cell_tot$proptot, method="spearman")

colors1=c("#3969AC","#80BA5A","#008695","#F2B701","#E73F74","#E68310","#CF1C90","#11A579","#f97b72","#A5AA99")

cell_tot$`Cell type`=gsub("_", " ", cell_tot$cell)
cell_tot$`Cell type`[which(cell_tot$`Cell type`=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"
cell_tot$`Cell type`=factor(cell_tot$`Cell type`, levels=c("Type 1", "Type 2a","Type 2x","Endothelial","Fibro-adipogenic progenitor","Smooth Muscle","T cell","Neuronal","Neuromuscular junction","Satellite Cell","Adipocyte","Macrophage"))

tiff("~/plot.tiff", height=180, width=250, units="mm", res=300)
ggplot(cell_tot, aes(x=total, y=proptot, color=`Cell type`))+geom_point(size=3)+theme_bw()+scale_color_manual(values=colors1)+theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=14))+scale_x_log10()+xlab("Total count")+ylab("Proportion genes differentially\nexpressed by sex (FDR<5%)")
dev.off()



