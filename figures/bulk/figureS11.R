library(data.table)
library(ggplot2)
library(gridExtra)


bulkcpm=fread("../snrna/cpms/bulk_cpm.tab")
pbulkcpm=fread("../snrna/cpms/pbulk_cpm.tab")
t2acpm=fread("../snrna/cpms/t2a_cpm.tab")
t2xcpm=fread("../snrna/cpms/t2x_cpm.tab")
bulkcpm$`Cell type`=rep("Bulk", nrow(bulkcpm))
pbulkcpm$`Cell type`=rep("Pseudobulk", nrow(pbulkcpm))
t2acpm$`Cell type`=rep("Type 2a", nrow(t2acpm))
t2xcpm$`Cell type`=rep("Type 2x", nrow(t2xcpm))
t2acol=grep("ENSG00000122705", colnames(t2acpm))
t2xcol=grep("ENSG00000122705", colnames(t2xcpm))
pcol=grep("ENSG00000122705", colnames(pbulkcpm))
fiber=rbind(as.data.frame(t2acpm)[,c(41438,t2acol,41439)], as.data.frame(t2xcpm)[,c(40528,t2xcol,40529)], as.data.frame(pbulkcpm)[,c(58963,pcol,58964)])
colnames(fiber)[2]="ENSG00000122705"
cpm=rbind(fiber, bulkcpm[,c("labelcode","ENSG00000122705","Cell type")])
pheno=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/tissue.csv")
cpm=merge(cpm, pheno[,c("labelcode","SEX")], by="labelcode")
print(max(cpm$ENSG00000122705))

cpm$pval=rep(1, nrow(cpm))
cpm$pval[which(cpm$`Cell type`=="Type 2a")]="p=0.005"
cpm$pval[which(cpm$`Cell type`=="Type 2x")]="p=0.003"
cpm$pval[which(cpm$`Cell type`=="Pseudobulk")]="p=0.001"
cpm$pval[which(cpm$`Cell type`=="Bulk")]="p=3x10-9"
cpm$facet=paste(cpm$`Cell type`, cpm$pval, sep="\n")

cpm$facet=factor(cpm$facet, levels=c("Type 2a\np=0.005", "Type 2x\np=0.003","Pseudobulk\np=0.001","Bulk\np=3x10-9"))

cpm=cpm[which(cpm$ENSG00000122705<200),]

tiff("~/plot.tiff", height=140, width=180, units="mm", res=300)
ggplot(cpm, aes(color=SEX,x=SEX, y=ENSG00000122705))+theme_bw()+geom_boxplot(fill=NA, outlier.shape=NA)+geom_jitter(size=0.05)+facet_wrap(.~facet, scales="free", nrow=1)+scale_color_manual(values=c("#e41a1c","#377eb8"))+theme(strip.background=element_rect(color=NA, fill=NA), axis.text.y=element_text(size=6, angle=45), axis.text.x=element_text(size=6),panel.spacing.x=unit(0,"in"),strip.text=element_text(size=7), legend.position="none", axis.title=element_blank(), plot.title=element_text(face="italic"), plot.margin=unit(c(0.07,0.07,0.07,0), "in"))
dev.off()




