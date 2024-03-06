library(data.table)


bulkpeak=as.data.frame(fread("/net/snowwhite/home/schanks/muscle/snatac/enrich/Bulk.1kb.peak.results"))
pbulkpeak=as.data.frame(fread("/net/snowwhite/home/schanks/muscle/snatac/enrich/Pbulk.1kb.peak.results"))
t1peak=as.data.frame(fread("/net/snowwhite/home/schanks/muscle/snatac/enrich/Type_1.1kb.peak.results"))
t2apeak=as.data.frame(fread("/net/snowwhite/home/schanks/muscle/snatac/enrich/Type_2a.1kb.peak.results"))
t2xpeak=as.data.frame(fread("/net/snowwhite/home/schanks/muscle/snatac/enrich/Type_2x.1kb.peak.results"))
bulkpeak$fdr=p.adjust(bulkpeak$pvalue, method="fdr")
pbulkpeak$fdr=p.adjust(pbulkpeak$pvalue, method="fdr")
t1peak$fdr=p.adjust(t1peak$pvalue, method="fdr")
t2apeak$fdr=p.adjust(t2apeak$pvalue, method="fdr")
t2xpeak$fdr=p.adjust(t2xpeak$pvalue, method="fdr")
bulkpeak$cell=rep("Bulk", nrow(bulkpeak))
pbulkpeak$cell=rep("Pseudobulk", nrow(pbulkpeak))


t1peak$male=(t1peak$fdr<0.05 & t1peak$log2FoldChange>0)
t1peak$female=(t1peak$fdr<0.05 & t1peak$log2FoldChange<0)
t2apeak$male=(t2apeak$fdr<0.05 & t2apeak$log2FoldChange>0)
t2apeak$female=(t2apeak$fdr<0.05 & t2apeak$log2FoldChange<0)
t2xpeak$male=(t2xpeak$fdr<0.05 & t2xpeak$log2FoldChange>0)
t2xpeak$female=(t2xpeak$fdr<0.05 & t2xpeak$log2FoldChange<0)
bulkpeak$male=(bulkpeak$fdr<0.05 & bulkpeak$log2FoldChange>0)
bulkpeak$female=(bulkpeak$fdr<0.05 & bulkpeak$log2FoldChange<0)
pbulkpeak$male=(pbulkpeak$fdr<0.05 & pbulkpeak$log2FoldChange>0)
pbulkpeak$female=(pbulkpeak$fdr<0.05 & pbulkpeak$log2FoldChange<0)

t1peak$M_peak_test=as.numeric(t1peak$M_sig_peak>0)
t1peak$M_peak_test[which(t1peak$F_sig_peak>0)]=NA
t1peak$F_peak_test=as.numeric(t1peak$F_sig_peak>0)
t1peak$F_peak_test[which(t1peak$M_sig_peak>0)]=NA
t2apeak$M_peak_test=as.numeric(t2apeak$M_sig_peak>0)
t2apeak$M_peak_test[which(t2apeak$F_sig_peak>0)]=NA
t2apeak$F_peak_test=as.numeric(t2apeak$F_sig_peak>0)
t2apeak$F_peak_test[which(t2apeak$M_sig_peak>0)]=NA
t2xpeak$M_peak_test=as.numeric(t2xpeak$M_sig_peak>0)
t2xpeak$M_peak_test[which(t2xpeak$F_sig_peak>0)]=NA
t2xpeak$F_peak_test=as.numeric(t2xpeak$F_sig_peak>0)
t2xpeak$F_peak_test[which(t2xpeak$M_sig_peak>0)]=NA
bulkpeak$M_peak_test=as.numeric(bulkpeak$M_sig_peak>0)
bulkpeak$M_peak_test[which(bulkpeak$F_sig_peak>0)]=NA
bulkpeak$F_peak_test=as.numeric(bulkpeak$F_sig_peak>0)
bulkpeak$F_peak_test[which(bulkpeak$M_sig_peak>0)]=NA
pbulkpeak$M_peak_test=as.numeric(pbulkpeak$M_sig_peak>0)
pbulkpeak$M_peak_test[which(pbulkpeak$F_sig_peak>0)]=NA
pbulkpeak$F_peak_test=as.numeric(pbulkpeak$F_sig_peak>0)
pbulkpeak$F_peak_test[which(pbulkpeak$M_sig_peak>0)]=NA

t1m=glm(t1peak$male~t1peak$M_peak_test, family="binomial")
summary(t1m)$coefficients
t2am=glm(t2apeak$male~t2apeak$M_peak_test, family="binomial")
summary(t2am)$coefficients
t2xm=glm(t2xpeak$male~t2xpeak$M_peak_test, family="binomial")
summary(t2xm)$coefficients
bulkm=glm(bulkpeak$male~bulkpeak$M_peak_test, family="binomial")
summary(bulkm)$coefficients
pbulkm=glm(pbulkpeak$male~pbulkpeak$M_peak_test, family="binomial")
summary(pbulkm)$coefficients

t1f=glm(t1peak$female~t1peak$F_peak_test, family="binomial")
summary(t1f)$coefficients
t2af=glm(t2apeak$female~t2apeak$F_peak_test, family="binomial")
summary(t2af)$coefficients
t2xf=glm(t2xpeak$female~t2xpeak$F_peak_test, family="binomial")
summary(t2xf)$coefficients
bulkf=glm(bulkpeak$female~bulkpeak$F_peak_test, family="binomial")
summary(bulkf)$coefficients
pbulkf=glm(pbulkpeak$female~pbulkpeak$F_peak_test, family="binomial")
summary(pbulkf)$coefficients



