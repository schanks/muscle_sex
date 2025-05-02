library(data.table)
library(ggplot2)

t1=fread("Type_1.1kb.peak.results")
t2a=fread("Type_2a.1kb.peak.results")
t2x=fread("Type_2x.1kb.peak.results")

t1$fdr=p.adjust(t1$pvalue, method="fdr")
t1=t1[which(t1$fdr<0.05),]
t1$Peaks=rep("0", nrow(t1))
t1$Peaks[which(t1$N_sig_peak>0)]="≥1"
t1$fc=abs(t1$log2FoldChange)
t1=t1[which(t1$chrom!="chrX" & t1$chrom!="chrY")]
t1$Peaks=factor(t1$Peaks, levels=c("0","≥1"))
t1$`Fiber type`=rep("Type 1", nrow(t1))

t2a$fdr=p.adjust(t2a$pvalue, method="fdr")
t2a=t2a[which(t2a$fdr<0.05),]
t2a$Peaks=rep("0", nrow(t2a))
t2a$Peaks[which(t2a$N_sig_peak>0)]="≥1"
t2a$fc=abs(t2a$log2FoldChange)
t2a=t2a[which(t2a$chrom!="chrX" & t2a$chrom!="chrY")]
t2a$Peaks=factor(t2a$Peaks, levels=c("0","≥1"))
t2a$`Fiber type`=rep("Type 2a", nrow(t2a))

t2x$fdr=p.adjust(t2x$pvalue, method="fdr")
t2x=t2x[which(t2x$fdr<0.05),]
t2x$Peaks=rep("0", nrow(t2x))
t2x$Peaks[which(t2x$N_sig_peak>0)]="≥1"
t2x$fc=abs(t2x$log2FoldChange)
t2x=t2x[which(t2x$chrom!="chrX" & t2x$chrom!="chrY")]
t2x$Peaks=factor(t2x$Peaks, levels=c("0","≥1"))
t2x$`Fiber type`=rep("Type 2x", nrow(t2x))

dat=rbind(t1, t2a, t2x)

tiff("~/plot.tiff", height=140, width=200, units="mm", res=300)
ggplot(dat, aes(x=Peaks, y=fc))+geom_boxplot(aes(fill=`Fiber type`))+theme_bw()+facet_grid(.~`Fiber type`)+theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="none", strip.text=element_text(size=14), strip.background=element_rect(fill=NA, color=NA))+scale_fill_manual(values=c("#3969AC","#80BA5A","#008695"))+ylab("Absolute value log2 fold change\nof autosomal sex-biased genes")+xlab("Number of sex-biased promoter peaks")
dev.off()

wilcox.test(t1$fc[which(t1$Peaks=="0")], t1$fc[which(t1$Peaks=="≥1")])
wilcox.test(t2a$fc[which(t2a$Peaks=="0")], t2a$fc[which(t2a$Peaks=="≥1")])
wilcox.test(t2x$fc[which(t2x$Peaks=="0")], t2x$fc[which(t2x$Peaks=="≥1")])



