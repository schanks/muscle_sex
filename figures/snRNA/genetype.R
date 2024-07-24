library(data.table)

celltypes=c("Type_1","Type_2a","Type_2x","Mesenchymal_Stem_Cell","Satellite_Cell","Neuromuscular_junction","Endothelial","Smooth_Muscle","Macrophage")

#Read in all results
for (c1 in celltypes){
	res=fread(paste("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/final_drop10nuc/results/",c1,".SEX.M.results.tab", sep=""))
	res$fdr=p.adjust(res$pvalue, method="fdr")
	if (c1=="Type_1"){
		dat=res
	}
	else{
		dat=rbind(dat, res)
	}
}

dat$sig=as.numeric(dat$fdr<0.05)
dat$cell=factor(dat$cell, levels=celltypes)
means=fread("cell_means_bygene.tab")
dat=merge(dat, means, by=c("gene", "cell"))
dat$type=rep("Other", nrow(dat))
dat$type[which(is.element(dat$gene_type, c("3prime_overlapping_ncRNA","antisense","bidirectional_promoter_lncRNA","lincRNA","macro_lncRNA","non_coding","processed_transcript","sense_intronic","sense_overlapping")))]="lncRNA"
dat$type[which(is.element(dat$gene_type, c("protein_coding")))]="Protein coding"
dat$type[which(is.element(dat$gene_type, c("polymorphic_pseudogene","pseudogene","processed_pseudogene","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","transcribed_unitary_pseudogene","unitary_pseudogene","unprocessed_pseudogene")))]="Pseudogene"
dat$type=factor(dat$type, levels=c("Protein coding","lncRNA","Pseudogene","Other"))
dat$meancat=cut(dat$meanCount, breaks=c(0,1,2,3,4,5,10,50,100,500,1000,5000,100000), labels=FALSE)
dat$meancat=as.factor(dat$meancat)

results=as.data.frame(matrix(c(rep(0, 6*9*3)), ncol=6))
colnames(results)=c("cell","genetype","OR","LB","UB","pvalue")
i=1
for (c1 in celltypes){
	fit=glm(sig~type+meancat,data=dat[which(dat$cell==c1),], family="binomial")
	res=summary(fit)$coefficients
	results$cell[i:(i+2)]=rep(c1,3)
	results$genetype[i:(i+2)]=rownames(res)[2:4]
	results$OR[i:(i+2)]=round(exp(as.numeric(res[2:4,1])), digits=2)
	results$LB[i:(i+2)]=round(exp(as.numeric(res[2:4,1])-1.96*(as.numeric(res[2:4,2]))),digits=2)
	results$UB[i:(i+2)]=round(exp(as.numeric(res[2:4,1])+1.96*(as.numeric(res[2:4,2]))),digits=2)
	results$pvalue[i:(i+2)]=signif(as.numeric(res[2:4,4]),2)
	i=i+3	
}


write.table(results,"genetype.tab", sep="\t", row.names=FALSE, quote=FALSE)



