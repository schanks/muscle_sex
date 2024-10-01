library(data.table)
library(ggplot2)
library(gridExtra)

#Panel A
noadj=fread("all.expr.results")
noadj=noadj[which(noadj$chrom!="chrY"),]
noadj$Chromosome=rep("Autosomal", nrow(noadj))
noadj$Chromosome[which(noadj$chrom=="chrX")]="X"
ogtt=fread("all.expr.ogtt.results")

celltypes=c("Type_1","Type_2a","Type_2x","Endothelial","Mesenchymal_Stem_Cell","Macrophage","Neuromuscular_junction","Neuronal","Satellite_Cell","Smooth_Muscle")

breaks=c(-100,-50,0,50,100)
labels=c("-100","-50","0","50","100")
sigcolors=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd")

for (c in celltypes){
	no=noadj[which(noadj$cell==c),c("gene","cell","Chromosome","log2FoldChange","pvalue","fdr")]
	og=ogtt[which(ogtt$cell==c),c("gene","log2FoldChange","pvalue","fdr")]
	both=merge(no, og, by="gene")
	colnames(both)=c("gene","cell","Chromosome","l2fc_no","p_no","fdr_no","l2fc_ogtt","p_ogtt","fdr_ogtt")
	both$log10p_no=-log10(both$p_no)
	both$log10p_no[which(both$l2fc_no<0)]=-both$log10p_no[which(both$l2fc_no<0)]
	both$log10p_ogtt=-log10(both$p_ogtt)
	both$log10p_ogtt[which(both$l2fc_ogtt<0)]=-both$log10p_ogtt[which(both$l2fc_ogtt<0)]
	both$Significance=rep("Neither", nrow(both))
	both$Significance[which(both$fdr_no<0.05 & both$fdr_ogtt<0.05)]="Both"
	both$Significance[which(both$fdr_no<0.05 & both$fdr_ogtt>0.05)]="No adjustment only"
	both$Significance[which(both$fdr_no>0.05 & both$fdr_ogtt<0.05)]="OGTT adjustment only"
	both$Significance=factor(both$Significance, levels=c("Both","No adjustment only","OGTT adjustment only","Neither"))
	bothsig=both[which(both$Significance=="Both")]
	concordance=(sum(bothsig$l2fc_no>0 & bothsig$l2fc_ogtt>0)+sum(bothsig$l2fc_no<0 & bothsig$l2fc_ogtt<0))/nrow(bothsig)
	concordance=round(concordance, digits=2)
	aut=bothsig[which(bothsig$Chromosome=="Autosomal"),]
	autc=(sum(aut$l2fc_no>0 & aut$l2fc_ogtt>0)+sum(aut$l2fc_no<0 & aut$l2fc_ogtt<0))/nrow(aut)
	autc=round(autc, digits=2)
	if (c=="Type_1"){ res=both}
	else {res=rbind(res, both)}
}

res$cell[which(res$cell=="Mesenchymal_Stem_Cell")]="Fibro-adipogenic progenitor"
res$cell=gsub("_"," ",res$cell)

a=ggplot(res, aes(x=log10p_no, y=log10p_ogtt))+facet_wrap(cell~., ncol=5)+ylab("Signed -log10 p-value: OGTT adjustment")+geom_point(aes(color=Significance, size=Chromosome))+theme_bw()+xlab("Signed -log10 p-value: No OGTT adjustment")+scale_y_continuous(breaks=breaks, labels=labels, limits=c(-1e+02, 1e+02))+scale_x_continuous(breaks=breaks, labels=labels, limits=c(-1e+02, 1e+02))+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+geom_abline(intercept=0, slope=1, linetype="dotted")+scale_shape_manual(values=c(16,17,18))+scale_size_manual(values=c(0.8,1.2,1.5))+scale_color_manual(values=sigcolors)+theme(legend.position="none", strip.background=element_rect(fill=NA, color=NA))

t1_noadj=fread("~/t1_out_invnorm.txt")
t1_ogtt=fread("~/t1_ogtt_out.txt")
t2a_noadj=fread("~/t2a_out_invnorm.txt")
t2a_ogtt=fread("~/t2a_ogtt_out.txt")
t2x_noadj=fread("~/t2x_out_invnorm.txt")
t2x_ogtt=fread("~/t2x_ogtt_out.txt")
t1=merge(t1_noadj[,c("Id","Name","#Genes","OddsRatio","P-Value","FDR")], t1_ogtt[,c("Id","OddsRatio","P-Value","FDR")], by="Id")
colnames(t1)=c("Id","Name","#Genes","OR_noadj","P_noadj","FDR_noadj","OR_ogtt","P_ogtt","FDR_ogtt")
t1$cell=rep("Type 1", nrow(t1))
t2a=merge(t2a_noadj[,c("Id","Name","#Genes","OddsRatio","P-Value","FDR")], t2a_ogtt[,c("Id","OddsRatio","P-Value","FDR")], by="Id")
colnames(t2a)=c("Id","Name","#Genes","OR_noadj","P_noadj","FDR_noadj","OR_ogtt","P_ogtt","FDR_ogtt")
t2a$cell=rep("Type 2a", nrow(t2a))
t2x=merge(t2x_noadj[,c("Id","Name","#Genes","OddsRatio","P-Value","FDR")], t2x_ogtt[,c("Id","OddsRatio","P-Value","FDR")], by="Id")
colnames(t2x)=c("Id","Name","#Genes","OR_noadj","P_noadj","FDR_noadj","OR_ogtt","P_ogtt","FDR_ogtt")
t2x$cell=rep("Type 2x",nrow(t2x))

all=rbind(t1, t2a, t2x)
all$log10p_noadj=-log10(all$P_noadj)
all$log10p_ogtt=-log10(all$P_ogtt)
all$log10p_noadj[which(all$OR_noadj<1)]=-all$log10p_noadj[which(all$OR_noadj<1)]
all$log10p_ogtt[which(all$OR_ogtt<1)]=-all$log10p_ogtt[which(all$OR_ogtt<1)]

all$Significance=rep("Neither", nrow(all))
all$Significance[which(all$FDR_noadj<0.05 & all$FDR_ogtt<0.05)]="Both"
all$Significance[which(all$FDR_noadj<0.05 & all$FDR_ogtt>0.05)]="No OGTT adjustment only"
all$Significance[which(all$FDR_noadj>0.05 & all$FDR_ogtt<0.05)]="OGTT adjustment only"
all$Significance=factor(all$Significance, levels=c("Both","No OGTT adjustment only","OGTT adjustment only","Neither"))

sigcols=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd")
b=ggplot(all, aes(x=log10p_noadj, y=log10p_ogtt, color=Significance))+geom_point(size=0.7)+scale_x_continuous(limits=c(-10.5,10.5))+scale_y_continuous(limits=c(-10.5,10.5))+theme_bw()+scale_color_manual(values=sigcols)+geom_abline(intercept=0, slope=1, linetype="dotted")+geom_hline(linetype="dotted",yintercept=0)+geom_vline(linetype="dotted",xintercept=0)+xlab("Signed -log10 p-value:\nNo OGTT adjustment")+ylab("Signed -log10 p-value:\nOGTT adjustment")+facet_grid(.~cell)+theme(legend.position="none",strip.background=element_rect(fill=NA, color=NA),plot.margin=unit(c(5.5,90,5.5,90),"pt"))

noadj=fread("../snatac/all.atac.results")
noadj$chrom=unlist(strsplit(noadj$peak, ":"))[seq(1, nrow(noadj)*3, by=3)]
noadj$Chromosome=rep("Autosomal", nrow(noadj))
noadj$Chromosome[which(noadj$chrom=="chrX")]="X"
ogtt=fread("../snatac/all.atac.ogtt.results")
celltypes=c("Type_1","Type_2a","Type_2x","Endothelial","Mesenchymal_Stem_Cell","Macrophage","Neuromuscular_junction","Neuronal","Satellite_Cell","Smooth_Muscle","T_cell","Adipocyte")

breaks=c(-100,-50,0,50,100)
labels=c("-100","-50","0","50","100")
sigcolors=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd")

for (c in celltypes){
	no=noadj[which(noadj$cell==c),c("peak","cell","Chromosome","log2FoldChange","pvalue","fdr")]
	og=ogtt[which(ogtt$cell==c),c("peak","log2FoldChange","pvalue","fdr")]
	both=merge(no, og, by="peak")
	colnames(both)=c("peak","cell","Chromosome","l2fc_no","p_no","fdr_no","l2fc_ogtt","p_ogtt","fdr_ogtt")
	both$log10p_no=-log10(both$p_no)
	both$log10p_no[which(both$l2fc_no<0)]=-both$log10p_no[which(both$l2fc_no<0)]
	both$log10p_ogtt=-log10(both$p_ogtt)
	both$log10p_ogtt[which(both$l2fc_ogtt<0)]=-both$log10p_ogtt[which(both$l2fc_ogtt<0)]
	both$Significance=rep("Neither", nrow(both))
	both$Significance[which(both$fdr_no<0.05 & both$fdr_ogtt<0.05)]="Both"
	both$Significance[which(both$fdr_no<0.05 & both$fdr_ogtt>0.05)]="No adjustment only"
	both$Significance[which(both$fdr_no>0.05 & both$fdr_ogtt<0.05)]="OGTT adjustment only"
	both$Significance=factor(both$Significance, levels=c("Both","No adjustment only","OGTT adjustment only","Neither"))
	bothsig=both[which(both$Significance=="Both")]
	concordance=(sum(bothsig$l2fc_no>0 & bothsig$l2fc_ogtt>0)+sum(bothsig$l2fc_no<0 & bothsig$l2fc_ogtt<0))/nrow(bothsig)
	concordance=round(concordance, digits=2)
	aut=bothsig[which(bothsig$Chromosome=="Autosomal"),]
	autc=(sum(aut$l2fc_no>0 & aut$l2fc_ogtt>0)+sum(aut$l2fc_no<0 & aut$l2fc_ogtt<0))/nrow(aut)
	autc=round(autc, digits=2)
	if (c=="Type_1"){ res=both}
	else {res=rbind(res, both)}
}

res$cell[which(res$cell=="Mesenchymal_Stem_Cell")]="Fibro-adipogenic progenitor"
res$cell=gsub("_"," ", res$cell)

c=ggplot(res, aes(x=log10p_no, y=log10p_ogtt))+facet_wrap(cell~., ncol=6)+ylab("Signed -log10 p-value: OGTT adjustment")+geom_point(aes(color=Significance))+theme_bw()+xlab("Signed -log10 p-value: No OGTT adjustment")+scale_y_continuous(breaks=breaks, labels=labels, limits=c(-1e+02, 1e+02))+scale_x_continuous(breaks=breaks, labels=labels, limits=c(-1e+02, 1e+02))+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+geom_abline(intercept=0, slope=1, linetype="dotted")+scale_shape_manual(values=c(16,17,18))+scale_size_manual(values=c(0.8,1.2,1.5))+scale_color_manual(values=sigcolors)+theme(legend.position="bottom", strip.background=element_rect(fill=NA, color=NA))+guides(color=guide_legend(override.aes=list(size=2)))

tiff("~/plot.tiff", units="in", width=11, height=13, res=200)
grid.arrange(a,b,c, heights=c(2,1.3,2))
dev.off()



