library(data.table)
library(ggplot2)
library(gridExtra)
library(scales)



colors1=c("#3969AC","#80BA5A","#008695","#F2B701","#E73F74","#E68310","#7F3C8D","#CF1C90","#11A579","#f97b72","#4b4b8f","#A5AA99")
colorsmfm=c("#3969AC","#80BA5A","#008695","#F2B701","#E73F74","#E68310","#7F3C8D","#CF1C90","#11A579","#f97b72","#4b4b8f","#A5AA99","#8c510a")

#Panel A: UMAP
nuc=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/hg38/fusion_plus_multiome_cluster_info/cluster_info_qc.tsv")
nuc=nuc[which(nuc$cohort=="FUSION"),]
nuc=nuc[which(nuc$modality=="atac"),]
exclude=omittedsampleids
nuc=nuc[which(!is.element(nuc$sample, exclude)),]
colnames(nuc)[17]="Cell type"
nuc$`Cell type`=gsub("_"," ",nuc$`Cell type`)
nuc$`Cell type`[which(nuc$`Cell type`=="Muscle Fiber Mixed")]="Mixed Muscle Fiber"
nuc$`Cell type`=factor(nuc$`Cell type`, levels=c("Type 1", "Type 2a","Type 2x","Endothelial","Mesenchymal Stem Cell","Smooth Muscle","T cell","Neuronal","Neuromuscular junction","Satellite Cell","Adipocyte","Macrophage","Mixed Muscle Fiber"))
a=ggplot(nuc, aes(x=UMAP_1, y=UMAP_2, color=`Cell type`))+theme_bw()+geom_point(size=0.5)+xlab("UMAP 1")+ylab("UMAP 2")+scale_color_manual(values=colorsmfm)+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.position="none")

#Panel B: Individual variability
dat=readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/data/df.atac.nuc.Rds")
dat[is.na(dat)]=0
dat$all=apply(dat[,-c(1:2)], 1, sum)
reprop=function(x){return(x/dat$all)}
props=apply(dat[,3:14], 2, reprop)
props=as.data.frame(props)
props$labelcode=dat$labelcode
pheno=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/tissue.csv")
pheno=pheno[,c("labelcode","SEX")]
dat=merge(props, pheno, by="labelcode")
colnames(dat)=gsub("_"," ",colnames(dat))
dat=dat[order(dat$SEX, dat$`Type 1`),]
datmelt=melt(dat[,-14], id.vars="labelcode")
colnames(datmelt)=c("labelcode","Cell type", "Proportion")
datmelt$labelcode=factor(datmelt$labelcode, levels=dat$labelcode)
mfm=datmelt[which(datmelt$`Cell type`=="Type 1"),]
mfm$`Cell type`=rep("Mixed Muscle Fiber", nrow(mfm))
mfm$Proportion=rep(0, nrow(mfm))
datmelt=rbind(datmelt, mfm)
datmelt$`Cell type`=as.character(datmelt$`Cell type`)
datmelt$`Cell type`[which(datmelt$`Cell type`=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"
datmelt$`Cell type`=factor(datmelt$`Cell type`, levels=c("Type 1", "Type 2a","Type 2x","Endothelial","Fibro-adipogenic progenitor","Smooth Muscle","T cell","Neuronal","Neuromuscular junction","Satellite Cell","Adipocyte","Macrophage","Mixed Muscle Fiber"))
b=ggplot(datmelt, aes(x=labelcode, fill=`Cell type`, y=Proportion))+geom_bar(position="stack",stat="identity", width=1)+theme_bw()+xlab(" ")+theme(axis.text.x=element_text(color="white"),axis.ticks.x=element_blank(), axis.text.y=element_text(size=12), axis.title=element_text(size=14))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=colorsmfm)+geom_text(label="Female:\nn=118", size=4.5,aes(x=5,y=0.94), check_overlap=TRUE, hjust=0, color="white")+geom_text(label="Male:\nn=163", size=4.5,aes(x=125,y=0.94), color="white",check_overlap=TRUE, hjust=0)

ab=grid.arrange(a, b, widths=c(1,2))

#Panel C - Mean proportions
fmeans=as.data.frame(as.numeric(apply(dat[which(dat$SEX=="F"),2:13],2,mean)))
colnames(fmeans)=c("Mean proportion")
fmeans$SD=as.numeric(apply(dat[which(dat$SEX=="F"),2:13],2,sd))
fmeans$Sex=rep("Female", nrow(fmeans))
fmeans$`Cell type`=colnames(dat)[2:13]
mmeans=as.data.frame(as.numeric(apply(dat[which(dat$SEX=="M"),2:13],2,mean)))
colnames(mmeans)=c("Mean proportion")
mmeans$SD=as.numeric(apply(dat[which(dat$SEX=="M"),2:13],2,sd))
mmeans$Sex=rep("Male", nrow(mmeans))
mmeans$`Cell type`=colnames(dat)[2:13]
means=rbind(fmeans, mmeans)
means$`Cell type`[which(means$`Cell type`=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"
means$`Cell type`=factor(means$`Cell type`, levels=c("Type 1", "Type 2a","Type 2x","Endothelial","Fibro-adipogenic progenitor","Smooth Muscle","T cell","Neuronal","Neuromuscular junction","Satellite Cell","Adipocyte","Macrophage"))
means$`Mean proportion`[which(means$Sex=="Female")]=-means$`Mean proportion`[which(means$Sex=="Female")]
means$textlabel=as.character(round(abs(means$`Mean proportion`), digits=2))
means$textlabel[which(means$textlabel=="0.1")]="0.10"
means$textlabel[which(means$textlabel=="0.2")]="0.20"
means$textpos=as.numeric(means$textlabel)+means$SD+0.08
means$textpos[which(means$Sex=="Female")]=-means$textpos[which(means$Sex=="Female")]
breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
labels=c("0.6","0.4","0.2","0.0","0.2","0.4","0.6")
c=ggplot(means, aes(x=`Mean proportion`, y=`Cell type`, fill=`Cell type`))+geom_bar(stat="identity")+geom_errorbar(aes(xmin=`Mean proportion`-SD, xmax=`Mean proportion`+SD),width=.3,linewidth=0.4)+scale_x_continuous(breaks=breaks, labels=labels, limits=c(-0.65,0.65))+scale_y_discrete(limits=rev)+theme_bw()+scale_fill_manual(values=colors1)+geom_text(aes(label=textlabel, x=textpos))+geom_text(label="Females", aes(x=-0.57, y=0.9), check_overlap=TRUE, size=4.5,hjust=0)+geom_vline(xintercept=0, size=0.2)+geom_text(label="Males", aes(x=.4, y=0.9), check_overlap=TRUE,size=4.5, hjust=0)+theme(legend.position="none", axis.text=element_text(size=12), axis.title=element_text(size=14))

#Panel D - Negative binomial results
dat=readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/output/nuclei_nb/nb_nuclei_invnorm_cell_SEX_fdr_across_celltypes_ATACnuc.Rds")
dat$OR=exp(dat$estimate)
dat$UB=exp(dat$estimate+1.96*dat$std.error)
dat$LB=exp(dat$estimate-1.96*dat$std.error)
dat$sig=rep("None", nrow(dat))
dat$sig[which(dat$estimate<0 & dat$padj<0.05)]="F"
dat$sig[which(dat$estimate>0 & dat$padj<0.05)]="M"
dat=dat[which(dat$celltype!="Muscle_Fiber_Mixed"),]
dat$celltype=gsub("_"," ", dat$celltype)
dat$celltype=factor(dat$celltype, levels=c("Type 1", "Type 2a","Type 2x","Endothelial","Mesenchymal Stem Cell","Smooth Muscle","T cell","Neuronal","Neuromuscular junction","Satellite Cell","Macrophage","Adipocyte"))
d=ggplot(dat, aes(y=celltype, x=OR, color=sig))+geom_vline(xintercept=1, linetype="longdash")+geom_point()+geom_errorbar(aes(xmin=LB, xmax=UB))+scale_y_discrete(limits=rev)+theme_bw()+scale_color_manual(values=c("#e41a1c","#377eb8","black"))+theme(legend.position="none", axis.title.x=element_text(size=14), axis.text.x=element_text(size=12), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+ylab("Cell type")+xlab("Fold change of male to female nuclei counts")+scale_x_continuous(trans=log_trans())

cd=grid.arrange(c,d, nrow=1)

tiff("~/plot.tiff", units="in", height=8, width=12, res=200)
grid.arrange(ab, cd, nrow=2)
dev.off()


