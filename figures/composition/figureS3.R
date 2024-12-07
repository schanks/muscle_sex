library(data.table)
library(ggplot2)
library(scales)

noadj=readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/output/nuclei_nb/nb_nuclei_invnorm_cell_SEX_fdr_across_celltypes.Rds")
ogtt=readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/output/nuclei_nb/nb_nuclei_invnorm_cell_fdrpercelltype_ogtt_covar.Rds")
ogtt=ogtt[which(ogtt$trait=="SEX.M"),]

celltypelevels=c("Type 1", "Type 2a","Type 2x","Endothelial","Mesenchymal Stem Cell","Smooth Muscle","T cell","Neuronal","Neuromuscular junction","Satellite Cell","Macrophage","Adipocyte")

noadj$OR=exp(noadj$estimate)
noadj$UB=exp(noadj$estimate+1.96*noadj$std.error)
noadj$LB=exp(noadj$estimate-1.96*noadj$std.error)
noadj$sig=rep("None", nrow(noadj))
noadj$sig[which(noadj$estimate<0 & noadj$padj<0.05)]="F"
noadj$sig[which(noadj$estimate>0 & noadj$padj<0.05)]="M"
noadj=noadj[which(noadj$celltype!="Muscle_Fiber_Mixed"),]
noadj$celltype=gsub("_"," ", noadj$celltype)
noadj$celltype=factor(noadj$celltype, levels=rev(celltypelevels))
noadj$Adjustment=rep("No OGTT adjustment", nrow(noadj))

ogtt$OR=exp(ogtt$estimate)
ogtt$UB=exp(ogtt$estimate+1.96*ogtt$std.error)
ogtt$LB=exp(ogtt$estimate-1.96*ogtt$std.error)
ogtt$sig=rep("None", nrow(ogtt))
ogtt$sig[which(ogtt$estimate<0 & ogtt$padj<0.05)]="F"
ogtt$sig[which(ogtt$estimate>0 & ogtt$padj<0.05)]="M"
ogtt=ogtt[which(ogtt$celltype!="Muscle_Fiber_Mixed"),]
ogtt$celltype=gsub("_"," ", ogtt$celltype)
ogtt$celltype=factor(ogtt$celltype, levels=rev(celltypelevels))
ogtt$Adjustment=rep("OGTT adjustment", nrow(ogtt))

both=rbind(noadj, ogtt)
both=as.data.frame(both)
both$yval=as.numeric(both$celltype)
both$yval[which(both$Adjustment=="OGTT adjustment")]=both$yval[which(both$Adjustment=="OGTT adjustment")]-0.2
both$yval[which(both$Adjustment=="No OGTT adjustment")]=both$yval[which(both$Adjustment=="No OGTT adjustment")]+0.2

tiff("~/FigureS3.tiff", height=4, width=6, units="in", res=200)
ggplot(both, aes(y=yval, x=OR, color=sig))+geom_vline(xintercept=1, linetype="longdash")+geom_point(size=1)+geom_errorbar(aes(xmin=LB, xmax=UB))+scale_y_continuous(breaks=c(seq(1,12)), labels=rev(celltypelevels))+theme_bw()+geom_point(aes(x=LB-0.02, shape=Adjustment))+ylab("Cell type")+scale_shape_manual(values=c(NA,8))+scale_color_manual(values=c("#e41a1c","#377eb8","black"))+theme(axis.title=element_text(size=7), axis.text=element_text(size=7))+xlab("Fold change of male to female nuclei counts")+scale_x_continuous(trans=log_trans())+guides(color="none")
dev.off()

