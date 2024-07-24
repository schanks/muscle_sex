library(ggplot2)
library(data.table)
library(gridExtra)

adjprop=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/bulk_adjprop/hg38.filter22323.bulk.SEX.M.results.tab")
noadj=fread("noadj/noadj.results.tab")
adjprop$fdr=p.adjust(adjprop$pvalue, method="fdr")
noadj$fdr=p.adjust(noadj$pvalue, method="fdr")
adjprop$log10p=-log10(adjprop$pvalue)
adjprop$log10p[which(adjprop$log2FoldChange<0)]=-adjprop$log10p[which(adjprop$log2FoldChange<0)]
noadj$log10p=-log10(noadj$pvalue)
noadj$log10p[which(noadj$log2FoldChange<0)]=-noadj$log10p[which(noadj$log2FoldChange<0)]

both=merge(noadj[,c("gene","symbol","log10p","fdr")], adjprop[,c("gene","log10p","fdr")], by="gene",all=TRUE)
colnames(both)[3:6]=c("p_no","f_no","p_adj","f_adj")
both$Significance=rep("Neither", nrow(both))
both$Significance[which(both$f_no<0.05 & both$f_adj<0.05)]="Both"
both$Significance[which(both$f_no<0.05 & both$f_adj>0.05)]="Unadjusted only"
both$Significance[which(both$f_no>0.05 & both$f_adj<0.05)]="Adjusted only"
both$Significance=factor(both$Significance, levels=c("Both","Unadjusted only","Adjusted only","Neither"))

breaks=c(-100,-50,0,50,100)
labels=c("-100","-50","0","50","100")
sigcols=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd")

a=ggplot(both, aes(x=p_no, y=p_adj, color=Significance))+geom_point(size=0.3)+xlab("Signed -log10 p-value:\nNo adjustment")+ylab("\n\nSigned -log10 p-value:\nAdjustment by estimated\ncell-type proporitons")+theme_bw()+scale_x_continuous(breaks=breaks, labels=labels,limits=c(-1e+02,1e+02))+scale_y_continuous(labels=labels,breaks=breaks, limits=c(-1e+02,1e+02))+theme(legend.title=element_text(size=8),legend.text=element_text(size=7),axis.title=element_text(size=7),axis.text=element_text(size=7))+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+geom_abline(intercept=0, slope=1, linetype="dotted")+scale_color_manual(values=sigcols)+guides(color=guide_legend(nrow=4, title.position="top", override.aes=list(size=1.5)))

t1count=fread("../snrna/cpms/t1_cpm.tab")
t1col=c(ncol(t1count),grep("ENSG00000092054", colnames(t1count)))
t1count=t1count[,..t1col]
colnames(t1count)=c("labelcode","cpm")
t1count$cell=rep("Type 1", nrow(t1count))
t2acount=fread("../snrna/cpms/t2a_cpm.tab")
t2acol=c(ncol(t2acount),grep("ENSG00000092054", colnames(t2acount)))
t2acount=t2acount[,..t2acol]
colnames(t2acount)=c("labelcode","cpm")
t2acount$cell=rep("Type 2a", nrow(t2acount))
t2xcount=fread("../snrna/cpms/t2x_cpm.tab")
t2xcol=c(ncol(t2xcount),grep("ENSG00000092054", colnames(t2xcount)))
t2xcount=t2xcount[,..t2xcol]
colnames(t2xcount)=c("labelcode","cpm")
t2xcount$cell=rep("Type 2x", nrow(t2xcount))
bulkcpm=fread("../snrna/cpms/bulk_cpm.tab")
bulkcol=c(ncol(bulkcpm),grep("ENSG00000092054", colnames(bulkcpm)))
bulkcpm=bulkcpm[,..bulkcol]
colnames(bulkcpm)=c("labelcode","cpm")
bulkcpm$cell=rep("Bulk", nrow(bulkcpm))
cpm=rbind(t1count, t2acount, t2xcount, bulkcpm)
pheno=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/tissue.csv")
cpm=merge(cpm, pheno[,c("labelcode","SEX")], by="labelcode")
cpm=as.data.frame(cpm)

b1=ggplot(cpm[which(cpm$cell=="Bulk"),], aes(color=SEX, x=SEX, y=cpm))+theme_bw()+geom_boxplot(fill=NA, outlier.shape=NA)+geom_jitter(size=0.1,height=0)+facet_grid(.~cell)+scale_color_manual(values=c("#e41a1c","#377eb8"))+theme(axis.title.y=element_text(size=7),strip.background=element_rect(color=NA, fill=NA), axis.text=element_text(size=6), strip.text=element_text(size=7), legend.position="none", axis.title.x=element_blank())+ylab("CPMs")

b2=ggplot(cpm[which(cpm$cell!="Bulk"),], aes(color=SEX, x=SEX, y=cpm))+theme_bw()+geom_boxplot(fill=NA, outlier.shape=NA)+geom_jitter(size=0.1,height=0)+facet_grid(.~cell)+scale_color_manual(values=c("#e41a1c","#377eb8"))+theme(axis.title.y=element_text(size=7),strip.background=element_rect(color=NA, fill=NA), axis.text=element_text(size=6), strip.text=element_text(size=7), legend.position="none", axis.title.x=element_blank())+ylab("CPMs")

btitle=ggplot(cpm, aes(x=1,y=1))+geom_text(size=3, hjust=0,check_overlap=TRUE, label="MYH7", fontface="italic")+scale_x_continuous(limits=c(1,10))+theme_void()+theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),plot.margin=unit(c(5.5,5.5,-8,5.5),"pt"))

b=grid.arrange(btitle, b1, b2, widths=c(1,2.5), heights=c(2,9), layout_matrix=rbind(c(1,1),c(2,3)))

tiff("~/plot.tiff", height=160, width=140, units="mm", res=300)
grid.arrange(a, b, heights=c(1,1))
dev.off()




