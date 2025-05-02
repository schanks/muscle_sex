library(data.table)
library(ggplot2)
library(gridExtra)

celltypes=c("Endothelial","Macrophage","Mesenchymal_Stem_Cell","Neuromuscular_junction","Neuronal","Satellite_Cell","Smooth_Muscle","Type_1","Type_2a","Type_2x")
cellorder=gsub("_"," ", celltypes)

for (c1 in celltypes[9]){
res1=fread(paste("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/output/DESeq.ATAC/final_drop10nuc/results/",c1,".SEX.M.results.tab", sep=""))
res1$fdr=p.adjust(res1$pvalue, method="fdr")
res1$r1sig=as.numeric(res1$fdr<0.05)
res1$log10p=-log10(res1$pvalue)
res1$log10p[which(res1$log2FoldChange<0)]=-res1$log10p[which(res1$log2FoldChange<0)]
res1$chrom=unlist(strsplit(res1$peak, ":"))[seq(1, nrow(res1)*3, by=3)]
res1=res1[,c("peak","chrom","log10p","r1sig")]
colnames(res1)=c("peak","chrom","FC_1","r1sig")
res1$Chromosome=rep("Autosomal", nrow(res1))
res1$Chromosome[which(res1$chrom=="chrX")]="X"

res2=NULL

for (c2 in celltypes){
	if (c2!=c1){
		res_c2=fread(paste("/net/snowwhite/home/aujackso/snRNAsnATAC_paper1/output/DESeq.ATAC/final_drop10nuc/results/",c2,".SEX.M.results.tab", sep=""))
		res_c2$fdr=p.adjust(res_c2$pvalue, method="fdr")
		res_c2$r2sig=as.numeric(res_c2$fdr<0.05)
		res_c2$log10p=-log10(res_c2$pvalue)
		res_c2$log10p[which(res_c2$log2FoldChange<0)]=-res_c2$log10p[which(res_c2$log2FoldChange<0)]
		res_c2=res_c2[,c("peak","log10p","cell","r2sig")]
		colnames(res_c2)=c("peak","FC_2","cell","r2sig")
		both=merge(res1[,c("peak","FC_1","r1sig","Chromosome")], res_c2[,c("peak","FC_2","r2sig")], by="peak")
		both=both[which(both$r1sig==1 & both$r2sig==1),]
		concordance=(sum(both$FC_1>0 & both$FC_2>0)+sum(both$FC_1<0 & both$FC_2<0))/nrow(both)
		res_c2$Concordance=rep(round(concordance, digits=2),nrow(res_c2))
		aut=both[which(both$Chromosome!="X" & both$Chromosome!="Y"),]
		autc=(sum(aut$FC_1>0 & aut$FC_2>0)+sum(aut$FC_1<0 & aut$FC_2<0))/nrow(aut)
		res_c2$AutConcordance=rep(round(autc, digits=2), nrow(res_c2))
		res2=rbind(res2, res_c2)
	}
}
res2$cell=gsub("_"," ", res2$cell)
res2$cell=factor(res2$cell, levels=cellorder)
res=merge(res1, res2, by="peak")
res$Significance=rep("Neither", nrow(res))
res$Significance[which(res$r1sig+res$r2sig==2)]="Both"
res$Significance[which(res$r1sig==1 & res$r2sig==0)]=paste(gsub("_"," ",c1),"only",sep=" ")
res$Significance[which(res$r1sig==0 & res$r2sig==1)]="Other cell type only"
res$Significance=factor(res$Significance, levels=c("Both",paste(gsub("_"," ",c1),"only",sep=" "),"Other cell type only","Neither"))

res$Concordance=as.character(res$Concordance)
res$Concordance[which(res$Concordance==1)]="1.00"
res$Concordance[which(res$Concordance==0.90)]="0.90"
res$Concordance[which(res$Concordance==NaN)]=" "
res$AutConcordance=as.character(res$AutConcordance)
res$AutConcordance[which(res$AutConcordance==1)]="1.00"
res$AutConcordance[which(res$AutConcordance==0.8)]="0.80"
res$AutConcordance[which(res$AutConcordance==0.9)]="0.90"
res$AutConcordance[which(res$AutConcordance==NaN)]=" "
res$Concordance=paste("All:\n", res$Concordance, sep="")
res$AutConcordance=paste("Autosomal:\n", res$AutConcordance, sep="")
res$Concordance[which(res$Concordance=="All:\n ")]=" "
res$AutConcordance[which(res$AutConcordance=="Autosomal:\n ")]=" "

cellname=gsub("_"," ",c1)
breaks=c(-100,-50,0,50,100)
labels=c("-100","-50","0","50","100")

res$cell=as.character(res$cell)
res$cell[which(res$cell=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"


tiff("~/plot.tiff", units="mm", width=160, height=180, res=300)
print(ggplot(res, aes(x=FC_1, y=FC_2, shape=Chromosome))+facet_wrap(cell~.,strip.position="left",nrow=4)+ggtitle(gsub("_"," ",c1))+ylab("Signed -log10 p-value:")+geom_point(aes(color=Significance,size=Chromosome))+geom_text(aes(label=Concordance,y=-90,x=100),hjust=1,size=2,check_overlap=TRUE)+geom_text(aes(label=AutConcordance,y=-45,x=100), hjust=1, size=2,check_overlap=TRUE)+theme_bw()+xlab(paste("Signed -log10 p-value: ",cellname))+scale_x_continuous(breaks=breaks, labels=labels,limits=c(-100,1e+02))+scale_y_continuous(labels=labels,breaks=breaks, limits=c(-100,100))+theme(axis.title=element_text(size=7),plot.title=element_text(size=7),legend.title=element_text(size=7),strip.placement="outside",legend.position="bottom",axis.text=element_text(size=7),strip.text=element_text(size=7),legend.text=element_text(size=7),strip.background=element_rect(fill="white",color=NA))+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+geom_abline(intercept=0, slope=1, linetype="dotted")+scale_shape_manual(values=c(16,17,18))+scale_size_manual(values=c(0.8,1.2,1.5))+scale_color_manual(values=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd"))+guides(color=guide_legend(nrow=4, title.position="top"),shape=guide_legend(nrow=3, override.aes=list(size=2), title.position="top")))
dev.off()

}

