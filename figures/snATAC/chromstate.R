library(data.table)

dat=fread("all.atac.results")
dat$sig=as.numeric(dat$fdr<0.05)
male=fread("chromstate/roadmap/peak_states_male.tab")
female=fread("chromstate/roadmap/peak_states_female.tab")
colnames(male)[4]="maleState"
colnames(female)[4]="femaleState"

dat=merge(dat, male[,c("peak","maleState")], by="peak")
dat=merge(dat, female[,c("peak","femaleState")], by="peak")

dat$chrom=unlist(strsplit(dat$peak, ":"))[seq(1, nrow(dat)*3, by=3)]
dat=dat[which(dat$chrom!="chrX"),]
dat=dat[which(dat$cell=="Type_1" | dat$cell=="Type_2a" | dat$cell=="Type_2x"),]

count=dat[,c("peak","maleState","femaleState")]
count=unique(count)
table(count$maleState, count$femaleState)
count$consensus=as.numeric(count$maleState==count$femaleState)

dat=dat[which(dat$maleState==dat$femaleState),]
dat$State=factor(dat$maleState, levels=c("15_Quies","1_TssA","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het"))

means=fread("cell_means_bypeak.tab")
dat=merge(dat, means, by=c("peak", "cell"))
dat$meancat=cut(dat$meanCount, breaks=c(0,1,2,3,4,5,10,50,100,1000), labels=FALSE)
dat$meancat=as.factor(dat$meancat)

state=dat[,c("peak","State")]
state=unique(state)
write.table(state, "consensus_peaks.tab", quote=FALSE, sep="\t", row.names=FALSE)

results=as.data.frame(matrix(c(rep(0, 6*3*14)), ncol=6))
colnames(results)=c("cell","State","OR","LB","UB","pvalue")
i=1
for (c1 in c("Type_1","Type_2a","Type_2x")){
	fit=glm(sig~State+meancat,data=dat[which(dat$cell==c1),], family="binomial")
	res=summary(fit)$coefficients
	results$cell[i:(i+13)]=rep(c1,14)
	results$State[i:(i+13)]=rownames(res)[2:15]
	results$OR[i:(i+13)]=round(exp(as.numeric(res[2:15,1])), digits=2)
	results$LB[i:(i+13)]=round(exp(as.numeric(res[2:15,1])-1.96*(as.numeric(res[2:15,2]))),digits=2)
	results$UB[i:(i+13)]=round(exp(as.numeric(res[2:15,1])+1.96*(as.numeric(res[2:15,2]))),digits=2)
	results$pvalue[i:(i+13)]=signif(as.numeric(res[2:15,4]),2)
	i=i+14
}

write.table(results,"15state_consensus.tab", sep="\t", row.names=FALSE, quote=FALSE)







