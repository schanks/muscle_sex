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
results=fread("/net/snowwhite/home/schanks/muscle/snrna/all.expr.results")
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
cell_tot=cell_tot[order(cell_tot$sigtot, decreasing=TRUE),]
cell_tot$cell=gsub("_"," ", cell_tot$cell)
cell_tot$cell[which(cell_tot$cell=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"
cell$cell=gsub("_"," ",cell$cell)
cell$cell[which(cell$cell=="Mesenchymal Stem Cell")]="Fibro-adipogenic progenitor"
cell$cell=factor(cell$cell, levels=c("Type 1","Type 2a","Type 2x","Fibro-adipogenic progenitor","Satellite Cell","Neuromuscular junction","Endothelial","Smooth Muscle","Macrophage","Neuronal"))
cell$legend=rep("F higher: Autosomal", nrow(cell))
cell$legend[which(cell$CHR=="ChrX" & cell$Direction==0)]="F higher: ChrX"
cell$legend[which(cell$CHR=="Autosomal" & cell$Direction==1)]="M higher: Autosomal"
cell$legend[which(cell$CHR=="ChrX" & cell$Direction==1)]="M higher: ChrX"
cell$legend[which(cell$CHR=="ChrY" & cell$Direction==1)]="M higher: ChrY"
cell$legend=as.factor(cell$legend)

down=fread("/net/snowwhite/home/schanks/muscle/snrna/all.down.fix.results")
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
cell$textheight=cell$proptot+0.005
cell$textheight[which(cell$cell=="Fibro-adipogenic progenitor")]=0.04

a=ggplot(cell, aes(x=cell, y=prop))+geom_bar(stat="identity",aes(fill=legend))+geom_segment(aes(y=downsample, yend=downsample, x=segstart, xend=segend, color=higher), linetype="dotted")+geom_text(aes(y=textheight, label=sigtot), size=2.2,check_overlap=TRUE)+theme_bw()+xlab(" ")+ylab("Proportion genes differentially\nexpressed by sex (FDR<5%)")+theme(plot.margin=unit(c(25,5.5,5.5,5.5),"pt"),legend.title=element_blank(), legend.position=c(0.7,0.75), axis.title=element_text(size=7),axis.text.y=element_text(size=7), axis.text.x=element_text(size=7, angle=45, hjust=0.95), legend.text=element_text(size=7),legend.margin=margin(unit(c(0,0,10,0),"mm")))+scale_fill_manual(values=c("#e41a1c","#67000d","#377eb8","#045a8d","#023858"))+guides(fill=guide_legend(keywidth=0.6, keyheight=0.6),title=element_blank())+scale_color_manual(values=c("white","black"))+guides(color="none")+geom_segment(y=0.08,yend=0.08, x=5.9,xend=6.3, linetype="dotted",check_overlap=TRUE)


#Panel C
dat=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/final_drop10nuc/results/Type_1.SEX.M.results.tab")
dat$sig=as.numeric(p.adjust(dat$pvalue, method="fdr")<0.05)
means=fread("/net/snowwhite/home/schanks/muscle/snrna/cell_means_bygene.tab")
dat=merge(dat, means, by=c("gene", "cell"))
dat$type=rep("Other", nrow(dat))
dat$type[which(is.element(dat$gene_type, c("3prime_overlapping_ncRNA","antisense","bidirectional_promoter_lncRNA","lincRNA","macro_lncRNA","non_coding","processed_transcript","sense_intronic","sense_overlapping")))]="lncRNA"
dat$type[which(is.element(dat$gene_type, c("protein_coding")))]="Protein coding"
dat$type[which(is.element(dat$gene_type, c("pseudogene","processed_pseudogene","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","transcribed_unitary_pseudogene","unitary_pseudogene","unprocessed_pseudogene")))]="Pseudogene"
dat$type=factor(dat$type, levels=c("lncRNA","Protein coding","Pseudogene","Other"))

c1=ggplot(dat[which(dat$type!="Other")], aes(x=meanCount))+geom_histogram(bins=50, aes(fill=type))+geom_text(aes(x=1000,y=450, label=type), size=2, hjust=0,check_overlap=TRUE)+scale_x_log10()+theme_bw()+facet_grid(type~.)+scale_fill_manual(values=c("#4daf4a","#984ea3","#ff7f00"))+scale_y_continuous(n.breaks=3)+theme(plot.margin=unit(c(25,5.5,0,4), "pt"),legend.position="none", strip.background=element_blank(), strip.text=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.spacing.y=unit(0, "mm"),axis.text.y=element_text(size=7), axis.title.y=element_text(size=7),axis.title.x=element_blank())+ylab("Count\n")

c2=ggplot(dat[which(dat$type!="Other"),], aes(x=meanCount, y=sig, group=type, color=type))+geom_smooth(span=0.01, se=TRUE, color=NA)+geom_smooth(span=0.01, se=FALSE)+geom_hline(yintercept=0)+geom_hline(yintercept=1)+scale_x_log10(breaks=c(0.1,10,1000,100000), labels=c("0.1","10","1K","100K"))+scale_y_continuous(limits=c(-0.13,1.11), breaks=c(0,0.5,1))+theme_bw()+xlab("Mean UMI in Type 1 fiber\n\n")+ylab("Proportion of genes differentially\nexpressed by sex (FDR<5%)")+theme(plot.margin=unit(c(0,0.07,0.07,0.07),"in"),legend.position=c(0.3,0.7),legend.background=element_rect(fill=NA, color=NA),axis.text.y=element_text(size=7),axis.text.x=element_text(size=7, hjust=0.95), axis.title=element_text(size=7), legend.title=element_text(size=7), legend.text=element_text(size=7), legend.key=element_rect(fill=NA, color=NA))+guides(color=guide_legend(title="Transcript type", keyheight=0.5))+scale_color_manual(values=c("#4daf4a","#984ea3","#ff7f00"))

c=grid.arrange(c1, c2, heights=c(1.5,3))

#Panel D
t1=fread("~/t1_out.txt")
t2a=fread("~/t2a_out.txt")
t2x=fread("~/t2x_out.txt")
t1=t1[which(t1$`#Genes`>=10),]
t2a=t2a[which(t2a$`#Genes`>=10),]
t2x=t2x[which(t2x$`#Genes`>=10),]
t1$cell=rep("Type 1", nrow(t1))
t2a$cell=rep("Type 2a", nrow(t2a))
t2x$cell=rep("Type 2x", nrow(t2x))
revigo=fread("/net/snowwhite/home/schanks/muscle/snrna/goterm/revigo.txt")
comp=rbind(t1[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")], t2a[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")], t2x[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")])
comp$signedp=-log10(comp$`P-Value`)
comp$signedp[which(comp$Direction=="down")]=-comp$signedp[which(comp$Direction=="down")]
comp=comp[which(is.element(comp$Name, revigo$Name)),]
minfdr=aggregate(comp$FDR, by=list(comp$Name), FUN=min)
minfdr$ormean=aggregate(comp$OddsRatio, by=list(comp$Name), FUN=mean)$x
minfdr=minfdr[minfdr$x<0.005,]
minfdr=minfdr[order(minfdr$ormean),]
comp=comp[is.element(comp$Name, minfdr$Group.1),]
comp$fdr=rep(NA, nrow(comp))
comp$fdr[which(comp$FDR<0.05)]="FDR<5%"
comp$fdr[which(comp$FDR<0.005)]="FDR<0.5%"
comp$fdr=factor(comp$fdr)

minfdr$Group.1[which(minfdr$Group.1=="positive regulation of morphogenesis of an epithelium")]="positive regulation of\nmorphogenesis of an epithelium"
minfdr$Group.1[which(minfdr$Group.1=="proton-transporting two-sector ATPase complex")]="proton-transporting two-sector\nATPase complex"
minfdr$Group.1[which(minfdr$Group.1=="regulation of smoothened signaling pathway")]="regulation of smoothened\nsignaling pathway"
minfdr$Group.1[which(minfdr$Group.1=="cellular response to growth factor stimulus")]="cellular response to\ngrowth factor stimulus"
minfdr$Group.1[which(minfdr$Group.1=="ubiquitin-dependent protein catabolic process")]="ubiquitin-dependent\nprotein catabolic process"
minfdr$Group.1[which(minfdr$Group.1=="energy derivation by oxidation of organic compounds")]="energy derivation by oxidation\nof organic compounds"
minfdr$Group.1[which(minfdr$Group.1=="generation of precursor metabolites and energy")]="generation of precursor\nmetabolites and energy"
comp$Name[which(comp$Name=="positive regulation of morphogenesis of an epithelium")]="positive regulation of\nmorphogenesis of an epithelium"
comp$Name[which(comp$Name=="proton-transporting two-sector ATPase complex")]="proton-transporting two-sector\nATPase complex"
comp$Name[which(comp$Name=="regulation of smoothened signaling pathway")]="regulation of smoothened\nsignaling pathway"
comp$Name[which(comp$Name=="cellular response to growth factor stimulus")]="cellular response to\ngrowth factor stimulus"
comp$Name[which(comp$Name=="ubiquitin-dependent protein catabolic process")]="ubiquitin-dependent\nprotein catabolic process"
comp$Name[which(comp$Name=="energy derivation by oxidation of organic compounds")]="energy derivation by oxidation\nof organic compounds"
comp$Name[which(comp$Name=="generation of precursor metabolites and energy")]="generation of precursor\nmetabolites and energy"
comp$Name=factor(comp$Name, levels=minfdr$Group.1)


d=ggplot(comp, aes(x=Name,y=cell, fill=log10(OddsRatio)))+geom_tile()+scale_fill_gradient2(low="#e41a1c",mid="white",high="#377eb8", midpoint=0,breaks=c(-0.6,-0.3,0,0.3), labels=c("-0.6","-0.3","0","0.3"))+geom_point(aes(shape=fdr), size=1)+theme_classic()+theme(legend.title=element_text(size=7),legend.position="bottom",axis.text.x=element_text(size=6), axis.text.y=element_text(size=7),axis.title=element_text(size=7),legend.text=element_text(size=7),plot.margin=unit(c(0.07,0,0.07,0.07),"in"))+xlab("GO Term")+ylab("Fiber type")+scale_shape_manual(values=c(8,1), na.translate=FALSE)+coord_flip()+guides(shape=guide_legend(direction="vertical", title=NULL, keyheight=0.8), fill=guide_colourbar(direction="horizontal", keyheight=1, keywidth=5, title.position="top"))

dlegend=get_legend(d)
d=d+theme(legend.position="none")
dwlegend=grid.arrange(d, dlegend,heights=c(5,1))

#Panel E= oxidative
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
layer_interval=10.8
oxires$testy=oxires$log10p+5
oxires$testy[which(oxires$cell=="Type_1")]=oxires$testy[which(oxires$cell=="Type_1")]+layer_interval*2
oxires$testy[which(oxires$cell=="Type_2a")]=oxires$testy[which(oxires$cell=="Type_2a")]+layer_interval
oxires$testend=rep(5, nrow(oxires))
oxires$testend[which(oxires$cell=="Type_1")]=5+layer_interval*2
oxires$testend[which(oxires$cell=="Type_2a")]=5+layer_interval*1
oxires$cell=gsub("_"," ", oxires$cell)
oxires=oxires[order(oxires$gene),]
oxires=oxires[which(is.element(oxires$gene, unique(oxires$gene)[1:50])),]
oxires$testy[which(is.element(oxires$gene, unique(oxires$gene)[41:50]))]=oxires$testend[which(is.element(oxires$gene, unique(oxires$gene)[41:50]))]

e=ggplot(oxires, aes(x=gene, y=testy, color=log2FoldChange>0))+geom_segment(aes(yend=testend, xend=gene))+geom_segment(y=5,yend=5,x=0,xend=40,size=0.3,color="gray")+geom_segment(y=5+layer_interval,yend=5+layer_interval,x=0,xend=40,size=0.3,color="gray")+geom_segment(y=5+layer_interval*2,yend=5+layer_interval*2,x=0,xend=40,size=0.3,color="gray")+theme_bw()+geom_text(x=0.5,y=5,label="0",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=5+layer_interval,label="0",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=5+layer_interval*2,label="0",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=10,label="5",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=10+layer_interval,label="5",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=10+layer_interval*2,label="5",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+theme(plot.margin=unit(c(0,-0.25,-0.5,0),"in"),legend.position="none",axis.text=element_blank(), axis.ticks=element_blank(), strip.background=element_rect(fill=NA, color=NA), panel.grid=element_blank(),panel.border=element_blank(),axis.title.y=element_blank())+scale_x_discrete(breaks=c())+scale_y_continuous(limits=c(0, 38.6),expand=expansion(mult=c(0,-0.1)))+geom_segment(x=0.5,xend=0.5,y=38.6,yend=5, size=0.3,color="gray")+geom_text(y=10+layer_interval,x=0.5,label="-log10 p-value",size=2.5,color="black",angle=90,check_overlap=TRUE, vjust=-1.5)+scale_color_manual(values=c("#e41a1c","#377eb8"))+xlab(" ")+geom_text(x=37.5,y=layer_interval*2+5,angle=90, check_overlap=TRUE, color="black",vjust=-0.5,size=3,label="GO:0006119\noxidative phosphorylation")+geom_text(aes(x=50, y=testend+5, label=cell), hjust=1.6,size=2.5,check_overlap=TRUE, color="black")+coord_polar()


#Panel G= caveola
CC=readRDS("/net/dumbo/home/dciotlos/FUSION/goCC.rds")
cav=CC$`GO:0005901`
cav=as.data.frame(cav)
colnames(cav)="Entrez_ID"
cavres=merge(cav, results, by="Entrez_ID")
cavres$log10p=-log10(cavres$pvalue)
cavres$log10p[which(cavres$log10p=="Inf")]=max(cavres$log10p[which(cavres$log10p!="Inf")])
mean_cav=aggregate(cavres$log10p, by=list(cavres$gene), FUN=mean)
mean_cav=mean_cav[order(mean_cav$x, decreasing=TRUE),]
cavres$gene=factor(cavres$gene, levels=c(mean_cav$`Group.1`))
layer_interval=17.9
cavres$testy=cavres$log10p+5
cavres$testy[which(cavres$cell=="Type_1")]=cavres$testy[which(cavres$cell=="Type_1")]+layer_interval*2
cavres$testy[which(cavres$cell=="Type_2a")]=cavres$testy[which(cavres$cell=="Type_2a")]+layer_interval
cavres$testend=rep(5, nrow(cavres))
cavres$testend[which(cavres$cell=="Type_1")]=5+layer_interval*2
cavres$testend[which(cavres$cell=="Type_2a")]=5+layer_interval
cavres$cell=gsub("_"," ", cavres$cell)
cavres=cavres[order(cavres$gene),]
cavres=cavres[which(is.element(cavres$gene, unique(cavres$gene)[1:50])),]
cavres$testy[which(is.element(cavres$gene, unique(cavres$gene)[41:50]))]=cavres$testend[which(is.element(cavres$gene, unique(cavres$gene)[41:50]))]

g=ggplot(cavres, aes(x=gene, y=testy, color=log2FoldChange>0))+geom_segment(y=5,yend=5,x=0,xend=40,size=0.3,color="gray")+geom_segment(y=5+layer_interval,yend=5+layer_interval,x=0,xend=40,size=0.3,color="gray")+geom_segment(y=5+layer_interval*2,yend=5+layer_interval*2,x=0,xend=40,size=0.3,color="gray")+geom_segment(aes(yend=testend, xend=gene))+theme_bw()+geom_text(x=0.5,y=5,label="0",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=5+layer_interval,label="0",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=5+layer_interval*2,label="0",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=10,label="5",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=10+layer_interval,label="5",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=10+layer_interval*2,label="5",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=15,label="10",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=15+layer_interval,label="10",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+geom_text(x=0.5,y=15+layer_interval*2,label="10",size=2,color="black",check_overlap=TRUE,hjust=1,vjust=0)+theme(plot.margin=unit(c(0,-0.2,-0.5,0),"in"),legend.position="none",axis.text=element_blank(), axis.ticks=element_blank(), strip.background=element_rect(fill=NA, color=NA), panel.grid=element_blank(),panel.border=element_blank(),axis.title.y=element_blank())+scale_x_discrete(breaks=c())+scale_y_continuous(limits=c(0, 58.7),expand=expansion(mult=c(0,-0.1)))+geom_segment(x=0.5,xend=0.5,y=58.7,yend=5, size=0.3,color="gray")+geom_text(y=10+layer_interval,x=0.5,label="-log10 p-value",size=2.5,color="black",angle=90,check_overlap=TRUE, vjust=-1.5)+scale_color_manual(values=c("#e41a1c","#377eb8"))+xlab(" ")+geom_text(x=37.5,y=2*layer_interval+5,angle=90,check_overlap=TRUE,color="black",size=3,vjust=-0.5,label="GO:0005901\ncaveola")+geom_text(aes(x=50, y=testend+5, label=cell), hjust=1.6,size=2.5,check_overlap=TRUE, color="black")+coord_polar()

eg=grid.arrange(e,g,nrow=2, heights=c(1,1))

row1=grid.arrange(a,c, widths=c(1.3,1))
row2=grid.arrange(dwlegend,eg, widths=c(2,1))

tiff("~/plot.tiff", height=180, width=180, units="mm", res=300)
grid.arrange(row1,row2,nrow=2, heights=c(1.8,2))
dev.off()

