library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpattern)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Panel A
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
both=both[which(both$Chromosome!="Y"),]

sigcolors=c("#ff7f00","#984ea3","#1b9e77","#bdbdbd")
a=ggplot(both, aes(x=bulksp, y=pbulksp, color=Significance, shape=Chromosome))+geom_point(size=0.8)+geom_vline(xintercept=0, linetype="dotted")+geom_hline(yintercept=0, linetype="dotted")+scale_x_continuous(limits=c(-100,100))+scale_y_continuous(limits=c(-100,100))+theme_bw()+scale_color_manual(values=sigcolors)+scale_shape_manual(values=c(16,2,5))+xlab("Bulk signed -log10 p-value")+ylab("Pseudobulk signed -log10 p-value")+theme(axis.title=element_text(size=7), axis.text=element_text(size=7), legend.title=element_text(size=7), legend.text=element_text(size=7), legend.position="right", plot.margin=unit(c(2,0.5,0.5,0.5), "lines"))+guides(color=guide_legend(nrow=4, keyheight=0.1, override.aes=list(size=1), title.position="top"), shape=guide_legend(keyheight=0.1,title.position="top",nrow=3, override.aes=list(size=1), order=2))

#Panel C
t1=fread("~/t1_out.txt")
t2a=fread("~/t2a_out.txt")
t2x=fread("~/t2x_out.txt")
t1=t1[which(t1$`#Genes`>=10),]
t2a=t2a[which(t2a$`#Genes`>=10),]
t2x=t2x[which(t2x$`#Genes`>=10),]
t1$cell=rep("Type 1", nrow(t1))
t2a$cell=rep("Type 2A", nrow(t2a))
t2x$cell=rep("Type 2X", nrow(t2x))
bulk=fread("~/bulk_out.txt")
bulk$cell=rep("Bulk", nrow(bulk))
pbulk=fread("~/pseudo_out.txt")
pbulk$cell=rep("P.bulk", nrow(pbulk))
revigo=fread("~/muscle/bulk/goterm/revigo_bulk.txt")
comp=rbind(t1[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")], t2a[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")], t2x[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")],bulk[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")],pbulk[,c("Name","Direction","P-Value","FDR","cell","OddsRatio")])
comp$signedp=-log10(comp$`P-Value`)
comp$signedp[which(comp$Direction=="down")]=-comp$signedp[which(comp$Direction=="down")]
comp=comp[which(is.element(comp$Name, revigo$Name)),]
minfdr=aggregate(comp$FDR, by=list(comp$Name), FUN=min)
minfdr$ormean=aggregate(comp$OddsRatio, by=list(comp$Name), FUN=mean)$x
minfdr=minfdr[minfdr$x<0.00000000005,]
minfdr=minfdr[order(minfdr$ormean),]
minfdr$Group.1[which(minfdr$Group.1=="nucleoside triphosphate metabolic process")]="nucleoside triphosphate\nmetabolic process"
minfdr$Group.1[which(minfdr$Group.1=="generation of precursor metabolites and energy")]="generation of precursor\nmetabolites and energy"
minfdr$Group.1[which(minfdr$Group.1=="glycosyl compound metabolic process")]="glycosyl compound\nmetabolic process"
minfdr$Group.1[which(minfdr$Group.1=="regulation of response to external stimulus")]="regulation of response\nto external stimulus"
minfdr$Group.1[which(minfdr$Group.1=="negative regulation of developmental process")]="negative regulation of\ndevelopmental process"
comp$Name[which(comp$Name=="nucleoside triphosphate metabolic process")]="nucleoside triphosphate\nmetabolic process"
comp$Name[which(comp$Name=="generation of precursor metabolites and energy")]="generation of precursor\nmetabolites and energy"
comp$Name[which(comp$Name=="glycosyl compound metabolic process")]="gl
ycosyl compound\nmetabolic process"
comp$Name[which(comp$Name=="regulation of response to external stimulus")]="regulation of response\nto external stimulus"
comp$Name[which(comp$Name=="negative regulation of developmental process")]="negative regulation of\ndevelopmental process"

comp=comp[is.element(comp$Name, minfdr$Group.1),]
comp$fdr=rep(NA, nrow(comp))
comp$fdr[which(comp$FDR<0.05)]="FDR<5%"
comp$fdr[which(comp$FDR<0.005)]="FDR<0.5%"
comp$fdr=factor(comp$fdr)
comp$Name=factor(comp$Name, levels=minfdr$Group.1)

c=ggplot(comp, aes(x=Name,y=cell, fill=log10(OddsRatio)))+geom_tile()+scale_fill_gradient2(low="#e41a1c",mid="white",high="#377eb8", midpoint=0,breaks=c(-0.2,-0.1,0,0.1,0.2,0.3), labels=c("-0.2","-0.1","0","0.1","0.2","0.3"))+geom_point(aes(shape=fdr), size=0.8)+theme_classic()+theme(legend.position="bottom",legend.title=element_text(size=7),axis.text.x=element_text(size=6, angle=45, hjust=0.95), axis.text.y=element_text(size=7),legend.box="vertical",axis.title.x=element_blank(), axis.title.y=element_text(size=7),legend.text=element_text(size=7),plot.margin=unit(c(0.2,0,0,0.2),"in"), legend.spacing.y=unit(0,"cm"))+xlab("GO Term")+scale_shape_manual(values=c(8,1), na.translate=FALSE)+coord_flip()+guides(shape=guide_legend(direction="vertical", title=NULL, keyheight=0.1, nrow=1), fill=guide_colourbar(direction="horizontal",barheight=0.5))

clegend=get_legend(c)
c=c+theme(legend.position="none")
cwlegend=grid.arrange(c, clegend,heights=c(4,1))

row1=grid.arrange(a, cwlegend,widths=c(1.4,1))

#Panels boxplots
bulk=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/bulk_adjprop/hg38.filter22323.bulk.SEX.M.results.tab")
bulk$fdr=p.adjust(bulk$pvalue, method="fdr")
snrna=fread("~/muscle/snrna/all.expr.results")
snrna$sig=as.numeric(snrna$fdr<0.05)
snrna$gene=unlist(strsplit(snrna$gene,"[.]"))[seq(1, nrow(snrna)*2, by=2)]
pseudo=fread("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.RNA/pbulk_adjprop/pb.SEX.M.results.tab")
pseudo$gene=unlist(strsplit(pseudo$gene,"[.]"))[seq(1, nrow(pseudo)*2, by=2)]
pseudo$fdr=p.adjust(pseudo$pvalue, method="fdr")
pseudo$sig=as.numeric(pseudo$fdr<0.05)
fiber=snrna[which(is.element(snrna$cell, c("Type_1","Type_2a","Type_2x"))),]
other=snrna[which(!is.element(snrna$cell, c("Type_1","Type_2a","Type_2x"))),]
fiber=fiber[which(fiber$fdr<0.05),]
other=other[which(other$fdr<0.05),]
fiberagg=aggregate(fiber$fdr, by=list(fiber$gene), FUN=min)
colnames(fiberagg)=c("gene","ffdr")
fiberagg$fmin=aggregate(fiber$log2FoldChange, by=list(fiber$gene), FUN=min)$x
fiberagg$fmax=aggregate(fiber$log2FoldChange, by=list(fiber$gene), FUN=max)$x
fiberagg$cells=aggregate(fiber$sig, by=list(fiber$gene), FUN=sum)$x
otheragg=aggregate(other$fdr, by=list(other$gene), FUN=min)
colnames(otheragg)=c("gene","ofdr")
otheragg$omin=aggregate(other$log2FoldChange, by=list(other$gene), FUN=min)$x
otheragg$omax=aggregate(other$log2FoldChange, by=list(other$gene), FUN=max)$x
snagg=merge(fiberagg, otheragg, by="gene", all=TRUE)
snagg=merge(snagg, bulk[,c("gene","log2FoldChange","fdr")], by="gene",all=TRUE)
colnames(snagg)[9:10]=c("bulkdir","bulkfdr")
snagg=merge(snagg, pseudo[,c("gene","log2FoldChange","fdr")], by="gene",all=TRUE)
colnames(snagg)[11:12]=c("pseudodir","pseudofdr")
snagg$ffdr[which(is.na(snagg$ffdr))]=100
snagg$ofdr[which(is.na(snagg$ofdr))]=100
snagg$bulkfdr[which(is.na(snagg$bulkfdr))]=100
snagg$pseudofdr[which(is.na(snagg$pseudofdr))]=100

#Panel C
snagg$Significance=rep("Not sex-biased in fiber types or bulk", nrow(snagg))
snagg$Significance[which(snagg$ffdr<0.05 & snagg$bulkfdr>0.05)]="Sex-biased in fiber types only"
snagg$Significance[which(snagg$ffdr>0.05 & snagg$bulkfdr<0.05)]="Sex-biased in bulk only"
snagg$Significance[which(snagg$ffdr<0.05 & snagg$bulkfdr<0.05)]="Sex-biased in fiber types and bulk"
snagg$Sig2=snagg$Significance
snagg$Sig2[which(snagg$Significance=="Sex-biased in fiber types and bulk" & sign(snagg$bulkdir)==sign(snagg$fmin) & sign(snagg$bulkdir)==sign(snagg$fmax))]="Concordant"
snagg$Sig2[which(snagg$Sig2=="Sex-biased in fiber types and bulk")]="Discordant"
snagg$Sig2=factor(snagg$Sig2, levels=c("Not sex-biased in fiber types or bulk","Sex-biased in bulk only","Sex-biased in fiber types only","Discordant","Concordant"))
snagg$Significance=factor(snagg$Significance, levels=c("Not sex-biased in fiber types or bulk","Sex-biased in bulk only","Sex-biased in fiber types only","Sex-biased in fiber types and bulk"))
cplot=aggregate(snagg$gene, by=list(snagg$Significance, snagg$Sig2), FUN=length)
colnames(cplot)=c("Significance","Sig2","count")

fplot=ggplot(cplot, aes(x=1.2,y=count, fill=Sig2))+geom_bar(width=0.4,stat="identity")+coord_flip()+theme_void()+theme(plot.margin=unit(c(5.5,5.5,5.5,0),"pt"),panel.grid=element_blank(),legend.position="none",axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+scale_x_continuous(limits=c(0.8,1.9))+scale_fill_manual(values=c("#bdbdbd","#984ea3","#1b9e77","#DB6D00","#FF9124"))+annotate("text",x=1.5,y=0, hjust=0,label="Fiber & bulk\n2,665", size=2.2)+annotate("text",x=1.5, y=2665, hjust=0,label="Fiber only\n2,463",size=2.2)+annotate("text",x=1.5,y=5150,hjust=0,label="Bulk only\n6,205",size=2.2)+annotate("text",x=1.5,y=11333,hjust=0,label="Neither\n14,067", size=2.2)+annotate("text",x=1.3,y=100,label="Concordant\nn=2,452",size=2,hjust=0)+annotate("text",x=1.1,y=2200,label="Discordant\nn=213",size=2,hjust=1)+annotate("segment",x=1.1,xend=1.1,y=2220,yend=2500,arrow=arrow(length=unit(3,"pt"),type="closed"))

textdat=as.data.frame(matrix(c("Comparison group","Bulk","Fiber","Bulk","Fiber","Fiber","Bulk","N/A","N genes",170,170,12,12,101,154,193,"N concordant",150,150,9,3,62,123,"N/A","% concordant","88%","88%","75%","25%","61%","80%","N/A"),byrow=TRUE, ncol=8))

oagg=snagg[which(snagg$ofdr<0.05),]
oagg$Concordance=rep("N/A", nrow(oagg))
oagg$Comparison=rep("Bulk", nrow(oagg))
oagg$Concordance[which(oagg$Sig2=="Concordant" & sign(oagg$bulkdir)==sign(oagg$omin) & sign(oagg$bulkdir)==sign(oagg$omax))]="Concordant"
oagg$Concordance[which(oagg$Sig2=="Concordant" & oagg$Concordance!="Concordant")]="Discordant"
oagg$Concordance[which(oagg$Sig2=="Discordant" & sign(oagg$bulkdir)==sign(oagg$omin) & sign(oagg$bulkdir)==sign(oagg$omax))]="Concordant"
oagg$Concordance[which(oagg$Sig2=="Discordant" & oagg$Concordance!="Concordant")]="Discordant"
oagg$Concordance[which(oagg$Sig2=="Sex-biased in fiber types only" & sign(oagg$omin)==sign(oagg$omax) & sign(oagg$omin)==sign(oagg$fmin) & sign(oagg$omin)==sign(oagg$fmax))]="Concordant"
oagg$Comparison[which(oagg$Sig2=="Sex-biased in fiber types only")]="Fiber type(s)"
oagg$Concordance[which(oagg$Sig2=="Sex-biased in fiber types only" & oagg$Concordance!="Concordant")]="Discordant"
oagg$Concordance[which(oagg$Sig2=="Sex-biased in bulk only" & sign(oagg$bulkdir)==sign(oagg$omin) & sign(oagg$bulkdir)==sign(oagg$omax))]="Concordant"
oagg$Concordance[which(oagg$Sig2=="Sex-biased in bulk only" & oagg$Concordance!="Concordant")]="Discordant"

conc=oagg[which(oagg$Sig2=="Concordant"),]
conc$Comparison=rep("Fiber type(s)", nrow(conc))

disc=oagg[which(oagg$Sig2=="Discordant"),]
disc$Concordance=rep("Discordant", nrow(disc))
disc$Concordance[which(sign(disc$fmax)==sign(disc$fmin) & sign(disc$omin)==sign(disc$fmin) & sign(disc$omax)==sign(disc$omin))]="Concordant"
disc$Comparison=rep("Fiber type(s)", nrow(disc))

fonly=oagg[which(oagg$Sig2=="Sex-biased in fiber types only"),]
fonly$Concordance=rep("N/A", nrow(fonly))
fonly$Comparison=rep("Bulk", nrow(fonly))

bonly=oagg[which(oagg$Sig2=="Sex-biased in bulk only"),]
bonly$Concordance=rep("N/A", nrow(bonly))
bonly$Comparison=rep("Fiber type(s)", nrow(bonly))

neither=oagg[which(oagg$Sig2=="Not sex-biased in fiber types or bulk"),]
neither$Concordance=rep("N/A", nrow(neither))
neither$Comparison=rep("Fiber type(s)", nrow(neither))

oagg=rbind(oagg, conc,disc, fonly, bonly, neither)


oagg$xpos=rep(900, nrow(oagg))
oagg$xpos[which(oagg$Sig2=="Concordant" & oagg$Comparison=="Bulk")]=1400
oagg$xpos[which(oagg$Sig2=="Discordant")]=2200
oagg$xpos[which(oagg$Sig2=="Discordant" & oagg$Comparison=="Bulk")]=2700
oagg$xpos[which(oagg$Sig2=="Sex-biased in fiber types only")]=3650
oagg$xpos[which(oagg$Sig2=="Sex-biased in fiber types only" & oagg$Comparison=="Bulk")]=4150
oagg$xpos[which(oagg$Sig2=="Sex-biased in bulk only")]=8000
oagg$xpos[which(oagg$Sig2=="Sex-biased in bulk only" & oagg$Comparison=="Bulk")]=8500
oagg$xpos[which(oagg$Sig2=="Not sex-biased in fiber types or bulk")]=18100
oagg$xpos[which(oagg$Sig2=="Not sex-biased in fiber types or bulk" & oagg$Comparison=="Bulk")]=18600

datplot=aggregate(oagg$gene, by=list(oagg$xpos, oagg$Comparison,oagg$Concordance), FUN=length)
colnames(datplot)=c("xpos","Comparison","Direction","Count")
datplot$Group=datplot$Comparison
datplot$Group[which(datplot$Direction=="N/A")]="N/A"
datplot$Direction[which(datplot$Group=="N/A")]="Concordant"

oplot=ggplot(datplot, aes(x=xpos, y=Count))+geom_col_pattern(position="fill",aes(color=Comparison,fill=Comparison,alpha=Group,pattern_density=Direction), pattern_color="black",pattern_size=0.2,pattern_fill=NA,pattern_spacing=0.04,pattern="crosshatch")+theme_bw()+scale_y_continuous(expand=expansion(add=c(0,0.1)))+theme(axis.title.y=element_text(size=7), axis.text.y=element_text(size=7),legend.title=element_text(size=7), legend.text=element_text(size=7),panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank(), panel.border=element_blank(),axis.text.x=element_text(size=6), axis.title.x=element_blank(), axis.ticks.x=element_blank(),legend.position="bottom")+ylab("Proportion concordant/discordant\nsex-biased expression")+scale_x_continuous(limits=c(0,25400),breaks=c(-1000,1150,2450,3900,8250,18350),labels=c("N non-fiber\ngenes","170","12","101","154","193"))+scale_fill_manual(values=c("#AD71B5","#49B192"))+guides(fill=guide_legend(keyheight=1,keywidth=1))+scale_pattern_density_manual(values=c(0,0.6,0))+scale_color_manual(values=c("#AD71B5","#49B192"))+scale_alpha_manual(values=c(1,1,0))+guides(alpha="none")+guides(pattern_density=guide_legend(override.aes=c(fill="white",color="black")))

legenddat=as.data.frame(matrix(c("Fiber","Concordant","Significant","Fiber","Discordant","Significant","Fiber","NA","Not significant","Bulk","Concordant","Significant","Bulk","Discordant","Significant","Bulk","NA","Not significant"),ncol=3,byrow=TRUE))
colnames(legenddat)=c("Comparison group","Directional concordance\nof non-fiber cell type","Fiber/bulk\nSignificance")
legenddat$`Comparison group`=factor(legenddat$`Comparison group`, levels=c("Fiber","Bulk"))
legenddat$`Directional concordance\nof non-fiber cell type`=factor(legenddat$`Directional concordance\nof non-fiber cell type`, levels=c("NA","Discordant","Concordant"))

olegend=ggplot(legenddat, aes(x=`Comparison group`,y=`Directional concordance\nof non-fiber cell type`))+geom_tile_pattern(aes(fill=`Comparison group`, alpha=`Fiber/bulk\nSignificance`, pattern_density=`Directional concordance\nof non-fiber cell type`,color=`Comparison group`),pattern_color="black",pattern_size=0.2,pattern_fill=NA,pattern_spacing=0.04,pattern="crosshatch",width=0.8,height=0.8,size=0.7)+scale_fill_manual(values=c("#49B192","#AD71B5"))+scale_color_manual(values=c("#49B192","#AD71B5"))+scale_alpha_manual(values=c(0,1))+scale_pattern_density_manual(values=c(0,0.6,0))+theme_void()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),legend.position="none", plot.margin=unit(c(0,5.5,-5.5,50),"pt"))+scale_x_discrete(expand=c(2,3))+scale_y_discrete(expand=c(0.05,1.5))+annotate("text",y=3.5,vjust=0,x=1,size=2,label="                  Comparison group:\nFiber-type")+annotate("text",y=3.5,vjust=0,x=2,size=2, label="Bulk")+annotate("text",y=3.5,vjust=0.2,x=-0.25,size=2.2,label="Direction\n")+annotate("text",y=3.5,vjust=0,x=-1.5,size=2,label="Comparison group\nsig.")+annotate("text",y=3,x=-0.25,size=2,label="Concordant")+annotate("text",y=2,x=-0.25,size=2,label="Discordant")+annotate("text",y=1,x=-0.25,size=2,label="N/A")+annotate("text",y=3,x=-1.5,size=2,label="Significant")+annotate("text",y=2,x=-1.5,size=2,label="Significant")+annotate("text",y=1,x=-1.5,size=2,label="Not significant")


oplot=oplot+theme(legend.position="none")
oplot_legend=ggplot(datplot, aes(x=xpos,y=2))+theme_void()+scale_x_continuous(limits=c(0,20000))
o=grid.arrange(oplot, olegend, heights=c(5,3))


textf=ggplot(legenddat, aes(x=-1,y=1))+theme_void()+annotate("text", hjust=0,size=2.2,vjust=0.6,y=1,x=1,label="Sex-biased expression significance\nand directional concordance\nof fiber types with bulk\nfor 25,400 tested genes")+theme(plot.margin=unit(c(5.5,5.5,5.5,-100),"pt"))

texto=ggplot(legenddat, aes(x=-1,y=1))+theme_void()+annotate("text",hjust=0,size=2.2,y=1,x=1,vjust=-0.5,label="Directional concordance of\n630 non-fiber cell type\nsex-biased genes\nwith fiber types (green)\nand bulk (purple)")+theme(plot.margin=unit(c(5.5,5.5,5.5,-100),"pt"))

textpanelc=ggplot(legenddat, aes(x=1,y=1))+theme_void()+annotate("text",hjust=0.5, size=2.2, y=2,x=2,vjust=-0.6, label="Sig. in:\nN genes")+theme(plot.margin=unit(c(5.5,-5.5,5.5,5.5),"pt"))


fplot=grid.arrange(textpanelc, fplot, widths=c(1,10))
row2plot=grid.arrange(fplot, o, heights=c(1,1.4))


row2text=grid.arrange(textf, texto, heights=c(1,1.4))

row2=grid.arrange(row2plot, row2text, widths=c(4,1))

tiff("~/plot.tiff", height=180, width=180, units="mm", res=300)
grid.arrange(row1, row2, heights=c(1,1.3))
dev.off()




