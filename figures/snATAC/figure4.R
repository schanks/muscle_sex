library(data.table)
library(ggplot2)
library(gridExtra)

sig_count<-function(x){
	return(sum(x<0.05))
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Panel A
results=fread("all.atac.results")
results$chrom=unlist(strsplit(results$peak, ":"))[seq(1, nrow(results)*3, by=3)]
results$CHR=rep("Autosomal", nrow(results))
results$CHR[which(results$chrom=="chrX")]="ChrX"
results$CHR[which(results$chrom=="chrY")]="ChrY"
results$Direction=as.numeric(results$log2FoldChange>0)

cell_tot=aggregate(results$fdr, by=list(results$cell), FUN=length)
colnames(cell_tot)=c("cell","tested")
cell_tot$sigtot=aggregate(results$fdr, by=list(results$cell), FUN=sig_count)$x
cell_tot$proptot=cell_tot$sigtot/cell_tot$tested
cell_prop=aggregate(results$fdr, by=list(results$cell, results$CHR, results$Direction), FUN=sig_count)
colnames(cell_prop)=c("cell","CHR","Direction","sig")
cell=merge(cell_tot, cell_prop, by="cell")
cell$prop=cell$sig/cell$tested
cell_tot=cell_tot[order(cell_tot$prop, decreasing=TRUE),]
cell_tot$cell=gsub("_"," ", cell_tot$cell)
cell_tot$cell[which(cell_tot$cell=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"
cell$cell=gsub("_"," ",cell$cell)
cell$cell[which(cell$cell=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"
cell$cell=factor(cell$cell, levels=c("Type 1","Type 2a","Type 2x","Endothelial","Fibro-adipogenic progenitor","Smooth Muscle","T cell","Neuronal","Neuromuscular junction","Satellite Cell","Adipocyte","Macrophage"))
cell$legend=rep("F higher: Autosomal", nrow(cell))
cell$legend[which(cell$CHR=="ChrX" & cell$Direction==0)]="F higher: ChrX"
cell$legend[which(cell$CHR=="Autosomal" & cell$Direction==1)]="M higher: Autosomal"
cell$legend[which(cell$CHR=="ChrX" & cell$Direction==1)]="M higher: ChrX"
cell$legend[which(cell$CHR=="ChrY" & cell$Direction==1)]="M higher: ChrY"
cell$legend=as.factor(cell$legend)
cell$sigtot[which(cell$sigtot==62680)]="63K"
cell$sigtot[which(cell$sigtot==54154)]="55K"
cell$sigtot[which(cell$sigtot==38197)]="39K"

a=ggplot(cell, aes(x=cell, y=prop))+geom_bar(stat="identity",aes(fill=legend))+geom_text(aes(y=proptot+0.005, label=sigtot), size=2,check_overlap=TRUE)+scale_y_continuous(limits=c(0,.135))+theme_bw()+xlab(" ")+ylab("Proportion peaks differentially\naccessible by sex")+theme(plot.margin=unit(c(0.3,0.07,0.07,0.07),"in"),legend.title=element_blank(), legend.position=c(0.65,0.75), axis.title=element_text(size=7),axis.text.y=element_text(size=7), legend.text=element_text(size=7),axis.text.x=element_text(size=7, angle=45, hjust=0.95))+scale_fill_manual(values=c("#e41a1c","#67000d","#377eb8","#045a8d"))+guides(fill=guide_legend(keyheight=0.5, keywidth=0.5))

#Panel B
down=fread("all.down.results")
down=down[which(down$baseMean!=0),]
down$cell[which(down$cell=="Mesenchymal_Stem_Cell")]="Fibro-adipogenic progenitor"
downcell=aggregate(down$fdr, by=list(down$cell), FUN=length)
colnames(downcell)=c("cell","tested")
downcell$sigtot=aggregate(down$fdr, by=list(down$cell), FUN=sig_count)$x
downcell$downsample=downcell$sigtot/downcell$tested
downcell$cell=gsub("_"," ",downcell$cell)

cell=merge(cell, downcell[,c("cell","downsample")], by="cell", all.x=TRUE)
cell$downsample[which(cell$cell=="Type 1")]=cell$proptot[which(cell$cell=="Type 1")]
cell$higher=as.factor(as.numeric(cell$downsample>cell$proptot))
cell$segstart=as.numeric(cell$cell)
cell$segstart=cell$segstart-0.4
cell$segend=cell$segstart+0.8

a=ggplot(cell, aes(x=cell, y=prop))+geom_bar(stat="identity",aes(fill=legend))+geom_segment(aes(x=segstart, xend=segend, y=downsample, yend=downsample, color=higher), linetype="dotted")+geom_text(aes(y=proptot+0.005, label=sigtot), size=2,check_overlap=TRUE)+scale_y_continuous(limits=c(0,.135))+theme_bw()+xlab(" ")+ylab("Proportion peaks differentially\naccessible by sex")+theme(plot.margin=unit(c(0.3,0.07,0.07,0.07),"in"),legend.title=element_blank(), legend.position=c(0.65,0.75), axis.title=element_text(size=7),axis.text.y=element_text(size=7), legend.text=element_text(size=7),axis.text.x=element_text(size=7, angle=45, hjust=0.95))+scale_fill_manual(values=c("#e41a1c","#67000d","#377eb8","#045a8d"))+guides(fill=guide_legend(keyheight=0.5, keywidth=0.5))+scale_color_manual(values=c("white","black"))+guides(color="none")


#Panel C
dat=results[which(results$cell=="Type_1"),]
dat$sig=as.numeric(dat$fdr<0.05)
means=fread("cell_means_bypeak.tab")
dat=merge(dat, means, by=c("peak","cell"))
state=fread("consensus_peaks.tab")
dat=merge(dat, state, by="peak")
dat=unique(dat)
dat=dat[which(is.element(dat$State, c("15_Quies","1_TssA","4_Tx","7_Enh"))),]
dat$State[which(dat$State=="15_Quies")]="Quiescent"
dat$State[which(dat$State=="1_TssA")]="Active TSS"
dat$State[which(dat$State=="4_Tx")]="Transcription"
dat$State[which(dat$State=="7_Enh")]="Enhancer"

c1=ggplot(dat, aes(x=meanCount))+geom_histogram(bins=50, aes(fill=State))+scale_x_log10()+theme_bw()+facet_grid(State~., scales="free_y")+scale_y_continuous(n.breaks=3)+theme(plot.margin=unit(c(20,5.5,0,7), "pt"),legend.position="none", strip.background=element_blank(), strip.text=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_text(size=7), axis.title.y=element_text(size=7, color="white"),axis.title.x=element_blank())+ylab("N peaks")

c2=ggplot(dat, aes(x=meanCount, y=sig, group=State, color=State))+geom_hline(yintercept=1)+geom_hline(yintercept=0)+geom_smooth(span=0.01, color=NA)+geom_smooth(span=0.01,se=FALSE)+scale_x_log10()+scale_y_continuous(limits=c(-0.05,1))+theme_bw()+xlab("Mean peak count\nin Type 1 fiber\n\n\n\n\n")+ylab("Proportion of peaks differentially\naccessible by sex")+theme(plot.margin=unit(c(0,5.5,5.5,5.5),"pt"),legend.position=c(0.3,0.7),axis.text.y=element_text(size=7),axis.text.x=element_text(size=7, hjust=0.95), axis.title.x=element_text(size=7), axis.title.y=element_text(size=7),legend.title=element_text(size=7), legend.text=element_text(size=7), legend.key=element_rect(fill=NA, color=NA),legend.background=element_blank())+guides(color=guide_legend(keywidth=0.1,keyheight=0.1, default.unit="inch"))+guides(shape=guide_legend(keywidth=0.1,keyheight=0.1, default.unit="inch"))

c=grid.arrange(c1, c2, heights=c(1.5,3))
row1=grid.arrange(a,c, widths=c(1.5,1))

#Panel D
celltypes=c("Type_1")
states=c("15_Quies","1_TssA","4_Tx","7_Enh")
tfbs_all=NULL
for (k in celltypes){
        for (s in states){
        d=fread(paste("tfbs/",k,"/TFBS_enrich_invnorm_results_",s,".tab", sep=""))
        tfbs_all=rbind(tfbs_all,d)
}
}
min=tfbs_all[which(tfbs_all$State=="15_Quies"),]
min=min[order(min$pvalue),]
min=min[1:5,]
minfdr=min
min=tfbs_all[which(tfbs_all$State=="1_TssA"),]
min=min[order(min$pvalue),]
min=min[1:5,]
minfdr=rbind(minfdr, min)
min=tfbs_all[which(tfbs_all$State=="4_Tx"),]
min=min[order(min$pvalue),]
min=min[1:5,]
minfdr=rbind(minfdr, min)
min=tfbs_all[which(tfbs_all$State=="7_Enh"),]
min=min[order(min$pvalue),]
min=min[1:5,]
minfdr=rbind(minfdr, min)
tfbs_all=tfbs_all[is.element(tfbs_all$TF, minfdr$TF),]
minfdr=aggregate(tfbs_all$OR, by=list(tfbs_all$TF), FUN=mean)
minfdr=minfdr[order(minfdr$x),]
tfbs_all$FDR=rep(NA, nrow(tfbs_all))
tfbs_all$FDR[which(tfbs_all$fdr<0.05)]="FDR<0.05"
tfbs_all$FDR[which(tfbs_all$fdr<5e-10)]="FDR<5e-10"
tfbs_all$FDR[which(tfbs_all$fdr<5e-100)]="FDR<5e-100"
tfbs_all$FDR=factor(tfbs_all$FDR)
tfbs_all$State[which(tfbs_all$State=="15_Quies")]="Quiescent"
tfbs_all$State[which(tfbs_all$State=="1_TssA")]="Active TSS"
tfbs_all$State[which(tfbs_all$State=="4_Tx")]="Strong\ntranscription"
tfbs_all$State[which(tfbs_all$State=="7_Enh")]="Enhancer"
states=c("Quiescent","Active TSS","Strong\ntranscription","Enhancer")

minfdr$Group.1[which(minfdr$Group.1=="CTCF_known2")]="CTCF\n{CTCF_known2}"
minfdr$Group.1[which(minfdr$Group.1=="NR3C1_known18")]="NR3C1\n{NR3C1_known18}"
minfdr$Group.1[which(minfdr$Group.1=="NR3C1_known6")]="AR\n{NR3C1_known6}"
minfdr$Group.1[which(minfdr$Group.1=="NR3C1_known13")]="AR\n{NR3C1_known13}"
minfdr$Group.1[which(minfdr$Group.1=="ZNF35_1")]="ZNF35\n{ZNF35_1}"
minfdr$Group.1[which(minfdr$Group.1=="PITX2_1")]="PITX2\n{PITX2_1}"
minfdr$Group.1[which(minfdr$Group.1=="NR3C1_known1")]="PGR\n{NR3C1_known1}"
minfdr$Group.1[which(minfdr$Group.1=="NR3C1_known11")]="NR3C1\n{NR3C1_known11}"
minfdr$Group.1[which(minfdr$Group.1=="PAX4_1")]="PAX4\n{PAX4_1}"
minfdr$Group.1[which(minfdr$Group.1=="PLAG1_1")]="PLAG1\n{PLAG1_1}"
minfdr$Group.1[which(minfdr$Group.1=="SPDEF_1")]="SPDEF\n{SPDEF_1}"
minfdr$Group.1[which(minfdr$Group.1=="FEV_1")]="FEV\n{FEV_1}"
minfdr$Group.1[which(minfdr$Group.1=="ELF3_3")]="ELF3\n{ELF3_3}"
minfdr$Group.1[which(minfdr$Group.1=="ETS_known7")]="ELK1\n{ETS_known7}"
minfdr$Group.1[which(minfdr$Group.1=="OVOL2_1")]="OVOL2\n{OVOL2_1}"
tfbs_all$TF[which(tfbs_all$TF=="CTCF_known2")]="CTCF\n{CTCF_known2}"
tfbs_all$TF[which(tfbs_all$TF=="NR3C1_known18")]="NR3C1\n{NR3C1_known18}"
tfbs_all$TF[which(tfbs_all$TF=="NR3C1_known6")]="AR\n{NR3C1_known6}"
tfbs_all$TF[which(tfbs_all$TF=="NR3C1_known13")]="AR\n{NR3C1_known13}"
tfbs_all$TF[which(tfbs_all$TF=="ZNF35_1")]="ZNF35\n{ZNF35_1}"
tfbs_all$TF[which(tfbs_all$TF=="PITX2_1")]="PITX2\n{PITX2_1}"
tfbs_all$TF[which(tfbs_all$TF=="NR3C1_known1")]="PGR\n{NR3C1_known1}"
tfbs_all$TF[which(tfbs_all$TF=="NR3C1_known11")]="NR3C1\n{NR3C1_known11}"
tfbs_all$TF[which(tfbs_all$TF=="PAX4_1")]="PAX4\n{PAX4_1}"
tfbs_all$TF[which(tfbs_all$TF=="PLAG1_1")]="PLAG1\n{PLAG1_1}"
tfbs_all$TF[which(tfbs_all$TF=="SPDEF_1")]="SPDEF\n{SPDEF_1}"
tfbs_all$TF[which(tfbs_all$TF=="FEV_1")]="FEV\n{FEV_1}"
tfbs_all$TF[which(tfbs_all$TF=="ELF3_3")]="ELF3\n{ELF3_3}"
tfbs_all$TF[which(tfbs_all$TF=="ETS_known7")]="ELK1\n{ETS_known7}"
tfbs_all$TF[which(tfbs_all$TF=="OVOL2_1")]="OVOL2\n{OVOL2_1}"


tfbs_all$TF=factor(tfbs_all$TF, levels=minfdr$Group.1)
tfbs_all$State=factor(tfbs_all$State, levels=c("Enhancer","Quiescent","Strong\ntranscription","Active TSS"))
tfbs_all$OR[which(tfbs_all$fdr>0.05)]=1


d=ggplot(tfbs_all, aes(x=TF, y=State, fill=log10(OR)))+geom_tile()+scale_fill_gradient2(low="#e41a1c",mid="white",high="#377eb8", midpoint=0, breaks=c(-0.05,0.05,0.15))+geom_point(aes(shape=FDR),size=1)+theme_bw()+theme(axis.text.x=element_text(size=6,angle=45, hjust=0.95), axis.text.y=element_text(size=7),axis.title=element_text(size=7),plot.margin=unit(c(0.07,0,0.07,0.07),"in"), legend.text=element_text(size=7), legend.title=element_text(size=7))+xlab("Transcription Factor: Type 1 fiber\n\n")+ylab("Chromatin state")+scale_shape_manual(values=c(1,16,8), na.translate=FALSE)+guides(shape=guide_legend(direction="vertical", title=NULL, keyheight=0.5), fill=guide_colourbar(direction="horizontal", title.position="top"))


#Panel E
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
enrich=rbind(t1peak[,c("chrom","log2FoldChange","fdr","F_sig_peak","M_sig_peak","cell")],t2apeak[,c("chrom","log2FoldChange","fdr","F_sig_peak","M_sig_peak","cell")],t2xpeak[,c("chrom","log2FoldChange","fdr","F_sig_peak","M_sig_peak","cell")],pbulkpeak[,c("chrom","log2FoldChange","fdr","F_sig_peak","M_sig_peak","cell")],bulkpeak[,c("chrom","log2FoldChange","fdr","F_sig_peak","M_sig_peak","cell")])
enrich=enrich[which(enrich$F_sig_peak<1  | enrich$M_sig_peak<1),]
enrich$cell=gsub("_"," ",enrich$cell)
enrich$cell=factor(enrich$cell, levels=c("Type 1","Type 2a","Type 2x","Pseudobulk","Bulk"))
enrich$`Sex-biased expression`=rep("Male higher; Significant")
enrich$`Sex-biased expression`[which(enrich$fdr>0.05 & enrich$log2FoldChange>0)]="Male higher; Not significant"
enrich$`Sex-biased expression`[which(enrich$fdr<0.05 & enrich$log2FoldChange<0)]="Female higher; Significant"
enrich$`Sex-biased expression`[which(enrich$fdr>0.05 & enrich$log2FoldChange<0)]="Female higher; Not significant"
enrich$`Sex-biased expression`=factor(enrich$`Sex-biased expression`, levels=c("Male higher; Significant","Male higher; Not significant","Female higher; Not significant","Female higher; Significant"))
enrich$Peakdir=enrich$M_sig_peak-enrich$F_sig_peak
enrich$Peakdir[which(enrich$Peakdir>1)]=1
enrich$Peakdir[which(enrich$Peakdir<(-1))]=-1
enrich$Peakdir=factor(enrich$Peakdir)
enrich=enrich[which(enrich$chrom!="chrX" & enrich$chrom!="chrY"),]

summary(as.factor(enrich$Peakdir))
sigcolors=c("#377eb8","#9BBFDC","#F28D8E","#e41a1c")

e1=ggplot(enrich, aes(fill=`Sex-biased expression`, y=Peakdir))+facet_grid(.~cell)+geom_bar(position="fill")+scale_fill_manual(values=sigcolors)+theme_bw()+xlab("Proportion of genes")+ylab("Number/direction of\nsex-biased promoter peaks")+scale_y_discrete(labels=c("≥1F","0","≥1M"))+scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.5","0.75","1"))+theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(), legend.position="bottom",strip.text=element_text(size=7),axis.text=element_text(size=7),axis.title=element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=7), plot.margin=unit(c(0.07,-0.02,0,0.07),"in"))+guides(fill=guide_legend(keyheight=0.7, keywidth=0.7, nrow=1, title.position="top"))

e2=ggplot(enrich[which(enrich$cell=="Bulk"),], aes(y=Peakdir))+facet_grid(.~cell)+geom_bar(stat="count")+theme_bw()+theme(axis.text.x=element_text(size=7), axis.title.x=element_text(size=7),panel.border=element_blank(),plot.margin=unit(c(0.07,0.07,0.07,0), "in"),strip.background=element_rect(fill=NA, color=NA), strip.text=element_text(color="white",size=7),axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.minor=element_blank(), panel.grid.major.y=element_blank())+scale_x_continuous(breaks=c(0,10000,20000), labels=c("0","10k","20k"))+xlab("N genes")

elegend=get_legend(e1)
e1=e1+theme(legend.position="none")
row3=grid.arrange(e1,e2,widths=c(5.5,1))

tiff("~/plot.tiff", height=180, width=180, units="mm", res=300)
grid.arrange(row1,d,row3,elegend,heights=c(2,1,0.8,0.25))
dev.off()





