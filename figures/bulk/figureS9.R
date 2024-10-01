library(data.table)
library(ggplot2)


bulk=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/bulk_adjprop/hg38.filter22323.bulk.SEX.M.results.tab")
gtex=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/gtex/hg38.filter22211.790gtex.SEX.M.results.tab")
gtex$gene=unlist(strsplit(gtex$gene,"[.]"))[seq(1, nrow(gtex)*2, by=2)]
bulk$ffdr=p.adjust(bulk$pvalue, method="fdr")
gtex$gfdr=p.adjust(gtex$pvalue, method="fdr")
bulk$bulk=-log10(bulk$pvalue)
bulk$bulk[which(bulk$log2FoldChange<0)]=-bulk$bulk[which(bulk$log2FoldChange<0)]
gtex$gtex=-log10(gtex$pvalue)
gtex$gtex[which(gtex$log2FoldChange<0)]=-gtex$gtex[which(gtex$log2FoldChange<0)]
bulk$Chromosome=rep("Autosomal",nrow(bulk))
bulk$Chromosome[which(bulk$chr=="X")]="X"
bulk=as.data.frame(bulk)
bulk=bulk[which(bulk$chr!="Y"),]

all=merge(bulk, gtex[,c("gene","gfdr","gtex")], by="gene")
all$Significance=rep("Neither", nrow(all))
all$Significance[which(all$ffdr<0.05 & all$gfdr<0.05)]="Both"
all$Significance[which(all$ffdr<0.05 & all$gfdr>0.05)]="FUSION only"
all$Significance[which(all$ffdr>0.05 & all$gfdr<0.05)]="GTEx only"
all$Significance=factor(all$Significance, levels=c("Both","FUSION only","GTEx only","Neither"))
both=all[which(all$Significance=="Both"),]
concordance=(sum(both$bulk>0 & both$gtex>0)+sum(both$bulk<0 & both$gtex<0))/nrow(both)
aut=both[which(both$Chromosome!="X" & both$Chromosome!="Y"),]
autc=(sum(aut$bulk>0 & aut$gtex>0)+sum(aut$bulk<0 & aut$gtex<0))/nrow(aut)
all$Concordance=rep(round(concordance, digits=2), nrow(all))
all$Concordance=paste("All:\n", all$Concordance, sep="")
all$AutConcordance=rep(round(autc, digits=2), nrow(all))
all$AutConcordance=paste("Autosomal:\n", all$AutConcordance, sep="")


breaks=c(-100,-50,0,50,100)
labels=c("-100","-50","0","50","100")

tiff("~/plot.tiff", units="in", width=4, height=5, res=300)
print(ggplot(all, aes(x=bulk, y=gtex, shape=Chromosome))+ylab("Signed -log10 p-value: GTEx")+xlab("Signed -log10 p-value: FUSION Bulk")+geom_point(aes(color=Significance,size=Chromosome))+geom_text(aes(label=Concordance,y=-90,x=100),hjust=1,check_overlap=TRUE)+geom_text(aes(label=AutConcordance, y=-45, x=100), hjust=1, check_overlap=TRUE)+theme_bw()+scale_x_continuous(breaks=breaks, labels=labels,limits=c(-1e+02,1e+02))+scale_y_continuous(labels=labels,breaks=breaks, limits=c(-1e+02,1e+02))+theme(plot.title=element_text(size=16),legend.title=element_text(size=10),strip.placement="outside",legend.position="bottom",axis.text=element_text(size=9),strip.text=element_text(size=10),strip.background=element_rect(fill="white",color=NA))+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+geom_abline(intercept=0, slope=1, linetype="dotted")+scale_shape_manual(values=c(16,17,18))+scale_size_manual(values=c(0.8,1.2,1.5))+scale_color_manual(values=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd"))+guides(color=guide_legend(nrow=4, title.position="top"),shape=guide_legend(nrow=3, override.aes=list(size=2), title.position="top")))
dev.off()

