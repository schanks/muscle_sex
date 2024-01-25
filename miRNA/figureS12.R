library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


bulk=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/bulk_adjprop/hg38.filter22323.bulk.SEX.M.results.tab")
bulk=bulk[which(bulk$baseMean!=0),]
bulk$bfdr=p.adjust(bulk$pvalue, method="fdr")
bulk$bulksp=-log10(bulk$pvalue)
bulk$bulksp[which(bulk$log2FoldChange<0)]=-bulk$bulksp[which(bulk$log2FoldChange<0)]
pseudo=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/pbulk_adjprop/pb.SEX.M.results.tab")
pseudo=pseudo[which(pseudo$baseMean!=0),]
pseudo$pfdr=p.adjust(pseudo$pvalue, method="fdr")
pseudo$pbulksp=-log10(pseudo$pvalue)
pseudo$pbulksp[which(pseudo$log2FoldChange<0)]=-pseudo$pbulksp[which(pseudo$log2FoldChange<0)]
pseudo$gene=unlist(strsplit(pseudo$gene,"[.]"))[seq(1, nrow(pseudo)*2, by=2)]
both=merge(bulk[,c("gene","bulksp","symbol","gene_type","chr","bfdr")], pseudo[,c("gene","pbulksp","pfdr")], by="gene")
both$Significance=rep("Neither", nrow(both))
both$Significance[which(both$bfdr<0.05 & both$pfdr>0.05)]="Bulk only"
both$Significance[which(both$bfdr>0.05 & both$pfdr<0.05)]="Pseudobulk only"
both$Significance[which(both$bfdr<0.05 & both$pfdr<0.05)]="Both"
both$Significance=factor(both$Significance, levels=c("Both","Bulk only","Pseudobulk only","Neither"))
both$Chromosome=rep("Autosomal", nrow(both))
both$Chromosome[which(both$chr=="X")]="X"
both$Chromosome[which(both$chr=="Y")]="Y"

#Panel F
mirna=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/miR/miR.filter.755.basecelltype.SEX.M.results.tab")
mirna$fdr=p.adjust(mirna$pvalue, method="fdr")
p5=mirna[grep("5p",mirna$gene),]
p3=mirna[grep("3p",mirna$gene),]
p5$transcript=unlist(strsplit(p5$gene,"5p"))
p3$transcript=unlist(strsplit(p3$gene,"3p"))
arms=merge(p5[,c("log2FoldChange","pvalue","fdr","chr","transcript")], p3[,c("log2FoldChange","pvalue","fdr","transcript")], by="transcript")
arms$Chromosome=rep("Autosomal",nrow(arms))
arms$Chromosome[which(arms$chr=="X")]="X"
arms$sp5=-log10(arms$pvalue.x)
arms$sp5[which(arms$log2FoldChange.x<0)]=-arms$sp5[which(arms$log2FoldChange.x<0)]
arms$sp3=-log10(arms$pvalue.y)
arms$sp3[which(arms$log2FoldChange.y<0)]=-arms$sp3[which(arms$log2FoldChange.y<0)]
arms$Significance=rep("Neither", nrow(arms))
arms$Significance[which(arms$fdr.x<0.05 & arms$fdr.y>0.05)]="5p only"
arms$Significance[which(arms$fdr.x>0.05 & arms$fdr.y<0.05)]="3p only"
arms$Significance[which(arms$fdr.x<0.05 & arms$fdr.y<0.05)]="Both"
arms$Significance=factor(arms$Significance, levels=c("Both","5p only","3p only","Neither"))
sigcolors=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd")

f=ggplot(arms, aes(x=sp5, y=sp3, color=Significance, shape=Chromosome))+geom_point(size=0.9)+geom_vline(xintercept=0, linetype="dotted")+geom_hline(yintercept=0, linetype="dotted")+scale_x_continuous(limits=c(-11,11))+scale_y_continuous(limits=c(-11,11))+theme_bw()+scale_color_manual(values=sigcolors)+xlab("5p arm signed -log10 p-value")+ylab("3p arm signed -log10 p-value")+theme(axis.title=element_text(size=7), axis.text=element_text(size=7), legend.title=element_text(size=7), legend.text=element_text(size=6), legend.position="bottom", plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))+guides(color=guide_legend(nrow=4, override.aes=list(size=1),title.position="top", keyheight=0.1), shape=guide_legend(nrow=4, keyheight=0.1, override.aes=list(size=1), title.position="top"))

#Panel G
target=read.csv('/net/snowwhite/home/asmauger/muscle/miRNAsexdifferences/mRNA_targets_df.csv', header=T, check.names=F, row.names = 1)
target$target=target$de_mir_targets_context60
target$target[which(target$target>=3)]=3
target$fdr=p.adjust(target$pvalue, method="fdr")
target$sig=as.numeric(target$fdr<0.05)
targetagg=aggregate(target$sig, by=list(target$target), FUN=mean)
colnames(targetagg)=c("target","prop")
targetagg$Genes=rep("All bulk", nrow(targetagg))
target=target[which(!is.element(target$gene, both$gene[which(both$pfdr<0.05)])),]
targetbulk=aggregate(target$sig, by=list(target$target), FUN=mean)
colnames(targetbulk)=c("target","prop")
targetbulk$Genes=rep("Bulk without\nsex-biased pseudobulk", nrow(targetbulk))
targetagg=rbind(targetagg, targetbulk)

g=ggplot(targetagg, aes(x=target, y=prop, color=Genes, group=Genes))+scale_color_manual(values=c("black","red"))+geom_point(size=0.8)+geom_line()+scale_y_continuous(limits=c(0,1))+theme_bw()+scale_x_continuous(breaks=c(0,1,2,3), labels=c("0","1","2","3+"))+xlab("Number of targeting\nsex-biased miRNAs")+ylab("Proportion of genes differentially\nexpressed by sex (FDR<5%)")+theme(axis.title=element_text(size=7), axis.text=element_text(size=7), panel.grid.minor.x=element_blank(), legend.text=element_text(size=7), plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"),legend.title=element_text(size=7), legend.position="bottom")+guides(color=guide_legend(nrow=2, title.position="top", keywidth=0.2, keyheight=0.8))


tiff("~/plot.tiff", height=80, width=100, units="mm", res=300)
grid.arrange(f, g, widths=c(1,1))
dev.off()




