library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

#Panel D - boxplot of a gene significant only in single nucleus data?
bulkcpm=fread("../snrna/cpms/bulk_cpm.tab")
pbulkcpm=fread("../snrna/cpms/pbulk_cpm.tab")
t2acpm=fread("../snrna/cpms/t2a_cpm.tab")
t2xcpm=fread("../snrna/cpms/t2x_cpm.tab")
bulkcpm$`Cell type`=rep("Bulk", nrow(bulkcpm))
pbulkcpm$`Cell type`=rep("Pseudobulk", nrow(pbulkcpm))
t2acpm$`Cell type`=rep("Type 2a", nrow(t2acpm))
t2xcpm$`Cell type`=rep("Type 2x", nrow(t2xcpm))
pheno=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/tissue.csv")

#Panel E - boxplot of a gene discordant from single nucleus to bulk
msccpm=fread("../snrna/cpms/msc_cpm.tab")
sccpm=fread("../snrna/cpms/sc_cpm.tab")
msccpm$`Cell type`=rep("Fibro-adipogenic\nprogenitor", nrow(msccpm))
sccpm$`Cell type`=rep("Satellite", nrow(sccpm))
msccol=grep("ENSG00000145012", colnames(msccpm))
sccol=grep("ENSG00000145012", colnames(sccpm))
t2acol=grep("ENSG00000145012", colnames(t2acpm))
t2xcol=grep("ENSG00000145012", colnames(t2xcpm))
pcol=grep("ENSG00000145012", colnames(pbulkcpm))
fiber=rbind(as.data.frame(msccpm)[,c(38001,msccol,38002)], as.data.frame(sccpm)[,c(32500,sccol,32501)], as.data.frame(pbulkcpm)[,c(58963,pcol,58964)], as.data.frame(t2acpm)[,c(41438,t2acol,41439)], as.data.frame(t2xcpm)[,c(40528,t2xcol,40529)])
colnames(fiber)[2]="ENSG00000145012"
cpm=rbind(fiber, bulkcpm[,c("labelcode","ENSG00000145012","Cell type")])
cpm=merge(cpm, pheno[,c("labelcode","SEX")], by="labelcode")

cpm$pval=rep(1, nrow(cpm))
cpm$pval[which(cpm$`Cell type`=="Type 2a")]="p=5x10-4"
cpm$pval[which(cpm$`Cell type`=="Type 2x")]="p=2x10-7"
cpm$pval[which(cpm$`Cell type`=="Fibro-adipogenic\nprogenitor")]="p=1x10-5"
cpm$pval[which(cpm$`Cell type`=="Satellite")]="p=3x10-6"
cpm$pval[which(cpm$`Cell type`=="Pseudobulk")]="0.15"
cpm$pval[which(cpm$`Cell type`=="Bulk")]="p=4x10-7"
cpm$facet=paste(cpm$`Cell type`, cpm$pval, sep="\n")

cpm$facet=factor(cpm$facet, levels=c("Type 2a\np=5x10-4", "Type 2x\np=2x10-7","Fibro-adipogenic\nprogenitor\np=1x10-5", "Satellite\np=3x10-6","Pseudobulk\n0.15","Bulk\np=4x10-7"))

tiff("~/plot.tiff", height=140, width=180, units="mm", res=300)
ggplot(cpm, aes(color=SEX,x=SEX, y=ENSG00000145012))+theme_bw()+geom_boxplot(fill=NA, outlier.shape=NA)+geom_jitter(size=0.05)+facet_wrap(.~facet, scales="free", nrow=1)+scale_color_manual(values=c("#e41a1c","#377eb8"))+theme(strip.background=element_rect(color=NA, fill=NA), axis.text.y=element_text(size=6, angle=45), axis.text.x=element_text(size=6),panel.spacing.x=unit(0,"in"),strip.text=element_text(size=7), legend.position="none", axis.title=element_blank(), plot.title=element_text(face="italic"), plot.margin=unit(c(0.07,0.07,0.07,0), "in"))
dev.off()






