library(data.table)

dat=fread("all.atac.results")
dat$sig=as.numeric(dat$fdr<0.05)
male=fread("chromstate/roadmap/peak_states_male.tab")

dat=merge(dat, male[,c("peak","State")], by="peak")

dat$chrom=unlist(strsplit(dat$peak, ":"))[seq(1, nrow(dat)*3, by=3)]
dat=dat[which(dat$chrom!="chrX"),]
dat=dat[which(dat$cell=="Type_1" | dat$cell=="Type_2a" | dat$cell=="Type_2x"),]

dat$State=factor(dat$State, levels=c("15_Quies","1_TssA","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het"))

results=as.data.frame(matrix(c(rep(0, 6*3*14)), ncol=6))
colnames(results)=c("cell","State","OR","LB","UB","pvalue")
i=1
for (c1 in c("Type_1","Type_2a","Type_2x")){
	fit=glm(sig~State,data=dat[which(dat$cell==c1),], family="binomial")
	res=summary(fit)$coefficients
	results$cell[i:(i+13)]=rep(c1,14)
	results$State[i:(i+13)]=rownames(res)[2:15]
	results$OR[i:(i+13)]=round(exp(as.numeric(res[2:15,1])), digits=2)
	results$LB[i:(i+13)]=round(exp(as.numeric(res[2:15,1])-1.96*(as.numeric(res[2:15,2]))),digits=2)
	results$UB[i:(i+13)]=round(exp(as.numeric(res[2:15,1])+1.96*(as.numeric(res[2:15,2]))),digits=2)
	results$pvalue[i:(i+13)]=signif(as.numeric(res[2:15,4]),2)
	i=i+14
}

write.table(results,"15state_nocount_male.tab", sep="\t", row.names=FALSE, quote=FALSE)






