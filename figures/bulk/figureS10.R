library(data.table)
library(ggplot2)

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

bulk_count=readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/data/bulk.hg38.counts.Rds")

pbulk_count=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/data/pbulk_counts.allgenes.txt")

meancpm<-function(x){
	x=as.numeric(x)
	return(mean(x))
}

bulkmean=as.data.frame(apply(bulk_count, 1, meancpm))
colnames(bulkmean)="bulkcpm"
bulk_count=as.data.frame(bulk_count)
bulkmean$gene=rownames(bulk_count)
pbulkmean=as.data.frame(apply(pbulk_count[,-1], 1, meancpm))
colnames(pbulkmean)="pbulkcpm"
pbulkmean$gene=unlist(strsplit(pbulk_count$gene,"[.]"))[seq(1, nrow(pbulk_count)*2, by=2)]

both=merge(both, pbulkmean, by="gene")
both=merge(both, bulkmean, by="gene")
both$pbulkbin=cut(both$pbulkcpm, breaks=c(-1,1,10,100,1000,10000,Inf), labels=FALSE)
both$bulkbin=cut(both$bulkcpm, breaks=c(-1,1,10,100,1000,10000,Inf), labels=FALSE)
both$pbulkbin[which(both$pbulkbin==1)]="Pseudobulk\ncount<1"
both$pbulkbin[which(both$pbulkbin==2)]="Pseudobulk\n1<count<10"
both$pbulkbin[which(both$pbulkbin==3)]="Pseudobulk\n10<count<100"
both$pbulkbin[which(both$pbulkbin==4)]="Pseudobulk\n100<count<1k"
both$pbulkbin[which(both$pbulkbin==5)]="Pseudobulk\n1k<count<10k"
both$pbulkbin[which(both$pbulkbin==6)]="Pseudobulk\ncount>10k"
both$bulkbin[which(both$bulkbin==1)]="Bulk\ncount<1"
both$bulkbin[which(both$bulkbin==2)]="Bulk\n1<count<10"
both$bulkbin[which(both$bulkbin==3)]="Bulk\n10<count<100"
both$bulkbin[which(both$bulkbin==4)]="Bulk\n100<count<1k"
both$bulkbin[which(both$bulkbin==5)]="Bulk\n1k<count<10k"
both$bulkbin[which(both$bulkbin==6)]="Bulk\ncount>10k"
both$pbulkbin=factor(both$pbulkbin, levels=c("Pseudobulk\ncount<1","Pseudobulk\n1<count<10","Pseudobulk\n10<count<100","Pseudobulk\n100<count<1k","Pseudobulk\n1k<count<10k","Pseudobulk\ncount>10k"))
both$bulkbin=factor(both$bulkbin, levels=c("Bulk\ncount<1","Bulk\n1<count<10","Bulk\n10<count<100","Bulk\n100<count<1k","Bulk\n1k<count<10k","Bulk\ncount>10k"))

tiff("~/plot.tiff", width=12, height=10, units="in", res=300)
ggplot(both[which(!is.element(both$chr, c("M","X","Y")))], aes(x=bulksp, y=pbulksp, color=Significance))+geom_point(size=0.8)+geom_vline(xintercept=0)+geom_hline(yintercept=0)+facet_grid(pbulkbin~bulkbin)+scale_x_continuous(limits=c(-10,10))+scale_y_continuous(limits=c(-10,10))+theme_bw()+scale_color_manual(values=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd"))+xlab("Bulk signed -log10 p-value")+ylab("Pseudobulk signed -log10 p-value")+theme(strip.background=element_rect(fill=NA, color=NA),strip.text=element_text(size=12), axis.title=element_text(size=14), axis.text=element_text(size=12), legend.title=element_text(size=14), legend.text=element_text(size=12))
dev.off()



