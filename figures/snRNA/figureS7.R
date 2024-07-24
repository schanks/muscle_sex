library(data.table)
library(ggplot2)
library(gridExtra)

results=fread("all.expr.results")
results$CHR=rep("Autosomal", nrow(results))
results$CHR[which(results$chrom=="chrX")]="ChrX"
results$CHR[which(results$chrom=="chrY")]="ChrY"
results$Direction=as.numeric(results$log2FoldChange>0)
BP=readRDS("/net/dumbo/home/dciotlos/FUSION/goBP.rds")
map=fread("/net/dumbo/home/dciotlos/R_Scripts/ensembl_to_entrez_V1_apr21.csv")
oxi=BP$`GO:0006119`
oxi=as.data.frame(oxi)
colnames(oxi)="Entrez_ID"
results=results[which(results$cell=="Type_1" | results$cell=="Type_2a" | results$cell=="Type_2x"),]
results$ENSEMBL=unlist(strsplit(results$gene, "[.]",))[seq(1, nrow(results)*2, by=2)]
results=merge(results, map[,c("ENSEMBL","Entrez_ID")], by="ENSEMBL")
results=results[which(results$chrom!="chrX" & results$chrom!="chrY"),]

oxires=merge(oxi, results, by="Entrez_ID")
oxires=unique(oxires)
oxires$log10p=-log10(oxires$pvalue)
oxires$log10p[which(oxires$log10p=="Inf")]=max(oxires$log10p[which(oxires$log10p!="Inf")])
mean_oxi=aggregate(oxires$log10p, by=list(oxires$gene), FUN=mean)
mean_oxi=mean_oxi[order(mean_oxi$x, decreasing=TRUE),]
oxires$gene=factor(oxires$gene, levels=c(mean_oxi$`Group.1`))
layer_interval=11.2
oxires$testy=oxires$log10p+5
oxires$testy[which(oxires$cell=="Type_1")]=oxires$testy[which(oxires$cell=="Type_1")]+layer_interval*2
oxires$testy[which(oxires$cell=="Type_2a")]=oxires$testy[which(oxires$cell=="Type_2a")]+layer_interval
oxires$testend=rep(5, nrow(oxires))
oxires$testend[which(oxires$cell=="Type_1")]=5+layer_interval*2
oxires$testend[which(oxires$cell=="Type_2a")]=5+layer_interval*1
oxires=oxires[order(oxires$gene),]
oxires=oxires[which(is.element(oxires$gene, unique(oxires$gene)[1:20])),]

genelist=unique(as.character(oxires$gene))

t1count=fread("cpms/t1_cpm.tab")
t1cols=c(which(is.element(colnames(t1count), genelist)),ncol(t1count))
t1count=t1count[,..t1cols]
t2acount=fread("cpms/t2a_cpm.tab")
t2acols=c(which(is.element(colnames(t2acount), genelist)),ncol(t2acount))
t2acount=t2acount[,..t2acols]
t2xcount=fread("cpms/t2x_cpm.tab")
t2xcols=c(which(is.element(colnames(t2xcount), genelist)),ncol(t2xcount))
t2xcount=t2xcount[,..t2xcols]
pheno=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/tissue.csv")
t1count=merge(t1count, pheno[,c("labelcode","SEX")], by="labelcode")
t2acount=merge(t2acount, pheno[,c("labelcode","SEX")], by="labelcode")
t2xcount=merge(t2xcount, pheno[,c("labelcode","SEX")], by="labelcode")
t1cpm=melt(t1count, id.vars=c("labelcode","SEX"))
t2acpm=melt(t2acount, id.vars=c("labelcode","SEX"))
t2xcpm=melt(t2xcount, id.vars=c("labelcode","SEX"))
t1cpm$`Fiber type`=rep("Type 1", nrow(t1cpm))
t2acpm$`Fiber type`=rep("Type 2a", nrow(t2acpm))
t2xcpm$`Fiber type`=rep("Type 2x", nrow(t2xcpm))
melted=rbind(t1cpm, t2acpm,t2xcpm)
melted$Sex=as.factor(melted$SEX)
melted$`Fiber type`=as.factor(melted$`Fiber type`)
colnames(melted)[3]="gene"


results$cell=gsub("_"," ",results$cell)
colnames(results)[9]="Fiber type"
melted=merge(melted, results[,c("gene","gene_name","Fiber type","pvalue")], by=c("gene","Fiber type"))
melted$pvalue=paste("p=",signif(melted$pvalue,2),sep="")

aplot=ggplot(melted, aes(x=`Fiber type`,y=value))+theme_bw()+geom_boxplot(aes(color=Sex),fill=NA, outlier.shape=NA,position=position_dodge(width=1))+scale_y_log10()+facet_wrap(gene_name~.,scales="free")+ylab("CPM")+scale_color_manual(values=c("#e41a1c","#377eb8"))+theme(panel.spacing.y=unit(0,"lines"),panel.spacing.x=unit(0.1, "lines"),strip.background=element_rect(color=NA, fill=NA), axis.text=element_text(size=7), strip.text.y=element_text(size=7),strip.text.x=element_text(size=7, face="italic"), axis.title=element_text(size=7), legend.position="none")
atitle=ggplot(melted, aes(x=1,y=1))+geom_text(size=3, hjust=0, check_overlap=TRUE, label="Oxidative phosphorylation")+scale_x_continuous(limits=c(1,10))+theme_void()+theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), plot.margin=unit(c(5.5,5.5,-8,5.5),"pt"))
a=grid.arrange(atitle, aplot, heights=c(1,11))

results=fread("all.expr.results")
results$CHR=rep("Autosomal", nrow(results))
results$CHR[which(results$chrom=="chrX")]="ChrX"
results$CHR[which(results$chrom=="chrY")]="ChrY"
results$Direction=as.numeric(results$log2FoldChange>0)
results=results[which(results$chrom!="chrX" & results$chrom!="chrY"),]
results=results[which(results$cell=="Type_1" | results$cell=="Type_2a" | results$cell=="Type_2x"),]
results$ENSEMBL=unlist(strsplit(results$gene, "[.]",))[seq(1, nrow(results)*2, by=2)]
CC=readRDS("/net/dumbo/home/dciotlos/FUSION/goCC.rds")

cav=CC$`GO:0005901`
cav=as.data.frame(cav)
colnames(cav)="Entrez_ID"
results=merge(results, map[,c("ENSEMBL","Entrez_ID")], by="ENSEMBL")
cavres=merge(cav, results, by="Entrez_ID")
cavres$log10p=-log10(cavres$pvalue)
cavres$log10p[which(cavres$log10p=="Inf")]=max(cavres$log10p[which(cavres$log10p!="Inf")])
mean_cav=aggregate(cavres$log10p, by=list(cavres$gene), FUN=mean)
mean_cav=mean_cav[order(mean_cav$x, decreasing=TRUE),]
cavres$gene=factor(cavres$gene, levels=c(mean_cav$`Group.1`))
cavres$cell=gsub("_"," ", cavres$cell)
cavres=cavres[order(cavres$gene),]
cavres=cavres[which(is.element(cavres$gene, unique(cavres$gene)[1:20])),]

genelist=unique(as.character(cavres$gene))

t1count=fread("cpms/t1_cpm.tab")
t1cols=c(which(is.element(colnames(t1count), genelist)),ncol(t1count))
t1count=t1count[,..t1cols]
t2acount=fread("cpms/t2a_cpm.tab")
t2acols=c(which(is.element(colnames(t2acount), genelist)),ncol(t2acount))
t2acount=t2acount[,..t2acols]
t2xcount=fread("cpms/t2x_cpm.tab")
t2xcols=c(which(is.element(colnames(t2xcount), genelist)),ncol(t2xcount))
t2xcount=t2xcount[,..t2xcols]
pheno=fread("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/data/tissue.csv")
t1count=merge(t1count, pheno[,c("labelcode","SEX")], by="labelcode")
t2acount=merge(t2acount, pheno[,c("labelcode","SEX")], by="labelcode")
t2xcount=merge(t2xcount, pheno[,c("labelcode","SEX")], by="labelcode")
t1cpm=melt(t1count, id.vars=c("labelcode","SEX"))
t2acpm=melt(t2acount, id.vars=c("labelcode","SEX"))
t2xcpm=melt(t2xcount, id.vars=c("labelcode","SEX"))
t1cpm$`Fiber type`=rep("Type 1", nrow(t1cpm))
t2acpm$`Fiber type`=rep("Type 2a", nrow(t2acpm))
t2xcpm$`Fiber type`=rep("Type 2x", nrow(t2xcpm))
melted=rbind(t1cpm, t2acpm,t2xcpm)
melted$Sex=as.factor(melted$SEX)
melted$`Fiber type`=as.factor(melted$`Fiber type`)
colnames(melted)[3]="gene"

results$cell=gsub("_"," ",results$cell)
colnames(results)[9]="Fiber type"
melted=merge(melted, results[,c("gene","gene_name","Fiber type","pvalue")], by=c("gene","Fiber type"))
melted$pvalue=paste("p=",signif(melted$pvalue,2),sep="")

bplot=ggplot(melted, aes(x=`Fiber type`,y=value))+theme_bw()+geom_boxplot(aes(color=Sex),fill=NA, outlier.shape=NA,position=position_dodge(width=1))+scale_y_log10()+facet_wrap(gene_name~.,scales="free")+ylab("CPM")+scale_color_manual(values=c("#e41a1c","#377eb8"))+theme(panel.spacing.y=unit(0,"lines"),panel.spacing.x=unit(0.1,"lines"),strip.background=element_rect(color=NA, fill=NA), axis.text=element_text(size=7), strip.text.y=element_text(size=7),strip.text.x=element_text(size=7, face="italic"), axis.title=element_text(size=7), legend.position="none")
btitle=ggplot(melted, aes(x=1,y=1))+geom_text(size=3, hjust=0, check_overlap=TRUE, label="Caveola")+scale_x_continuous(limits=c(1,10))+theme_void()+theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), plot.margin=unit(c(5.5,5.5,-8,5.5),"pt"))
b=grid.arrange(btitle, bplot, heights=c(1,11))

tiff("~/plot.tiff", units="mm", height=180, width=180, res=300)
grid.arrange(a,b, heights=c(1,1))
dev.off()

