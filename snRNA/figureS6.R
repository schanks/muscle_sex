library(data.table)
library(ggplot2)


dat=fread("genetype.tab")
dat$cell=gsub("_"," ",dat$cell)
dat$UB[which(dat$UB=="Inf")]=NA
dat$LB[which(dat$LB==0)]=NA
dat$genetype[which(dat$genetype=="typelncRNA")]="lncRNA"
dat$genetype[which(dat$genetype=="typeOther")]="Other"
dat$genetype[which(dat$genetype=="typePseudogene")]="Pseudogene"
dat$cell[which(dat$cell=="Mesenchymal Stem Cell")]="FAP"
dat$cell[which(dat$cell=="Neuromuscular junction")]="NMJ"
dat$cell[which(dat$cell=="Smooth Muscle")]="Smooth\nMuscle"
dat$genetype=factor(dat$genetype, levels=c("lncRNA","Pseudogene","Other"))
dat$Genes=rep("All genes", nrow(dat))

aut=fread("genetype_aut.tab")
aut$cell=gsub("_"," ",aut$cell)
aut$UB[which(aut$UB=="Inf")]=NA
aut$LB[which(aut$LB==0)]=NA
aut$genetype[which(aut$genetype=="typelncRNA")]="lncRNA"
aut$genetype[which(aut$genetype=="typeOther")]="Other"
aut$genetype[which(aut$genetype=="typePseudogene")]="Pseudogene"
aut$cell[which(aut$cell=="Mesenchymal Stem Cell")]="FAP"
aut$cell[which(aut$cell=="Neuromuscular junction")]="NMJ"
aut$cell[which(aut$cell=="Smooth Muscle")]="Smooth\nMuscle"
aut$genetype=factor(aut$genetype, levels=c("lncRNA","Pseudogene","Other"))
aut$Genes=rep("Autosomal genes", nrow(aut))

dat=rbind(aut, dat)
dat$OR[which(dat$OR==0)]=NA
dat$Adjustment=rep("UMI adjusted", nrow(dat))

d=fread("genetype_nocount.tab")
d$cell=gsub("_"," ",d$cell)
d$UB[which(d$UB=="Inf")]=NA
d$LB[which(d$LB==0)]=NA
d$genetype[which(d$genetype=="typelncRNA")]="lncRNA"
d$genetype[which(d$genetype=="typeOther")]="Other"
d$genetype[which(d$genetype=="typePseudogene")]="Pseudogene"
d$cell[which(d$cell=="Mesenchymal Stem Cell")]="FAP"
d$cell[which(d$cell=="Neuromuscular junction")]="NMJ"
d$cell[which(d$cell=="Smooth Muscle")]="Smooth\nMuscle"
d$genetype=factor(d$genetype, levels=c("lncRNA","Pseudogene","Other"))
d$Genes=rep("All genes", nrow(d))

a=fread("genetype_nocount_aut.tab")
a$cell=gsub("_"," ",a$cell)
a$UB[which(a$UB=="Inf")]=NA
a$LB[which(a$LB==0)]=NA
a$genetype[which(a$genetype=="typelncRNA")]="lncRNA"
a$genetype[which(a$genetype=="typeOther")]="Other"
a$genetype[which(a$genetype=="typePseudogene")]="Pseudogene"
a$cell[which(a$cell=="Mesenchymal Stem Cell")]="FAP"
a$cell[which(a$cell=="Neuromuscular junction")]="NMJ"
a$cell[which(a$cell=="Smooth Muscle")]="Smooth\nMuscle"
a$genetype=factor(a$genetype, levels=c("lncRNA","Pseudogene","Other"))
a$Genes=rep("Autosomal genes", nrow(a))

nocount=rbind(a, d)
nocount$OR[which(nocount$OR==0)]=NA
nocount$Adjustment=rep("No adjustment", nrow(nocount))

dat=rbind(dat, nocount)
dat$xval=rep(1, nrow(dat))
dat$xval[which(dat$genetype=="lncRNA")]=3
dat$xval[which(dat$genetype=="Pseudogene")]=2
dat$xval[which(dat$Adjustment=="No adjustment")]=dat$xval[which(dat$Adjustment=="No adjustment")]+0.3

tiff("~/plot.tiff", units="mm",height=100, width=120, res=300)
ggplot(dat, aes(y=xval))+theme_bw()+geom_vline(xintercept=1, linetype="longdash")+geom_segment(linewidth=0.8,aes(linetype=Adjustment, color=genetype,yend=xval,x=LB, xend=UB))+geom_point(size=1,aes(x=OR, color=genetype, shape=Adjustment), fill="white")+scale_x_log10(limits=c(0.03,160), position="top")+scale_shape_manual(values=c(19,21))+facet_grid(cell~Genes,switch="both")+theme(axis.title=element_text(size=7), axis.text.x=element_text(size=7),strip.text.y.left=element_text(size=7, angle=0),strip.text.x=element_text(size=7),legend.title=element_text(size=7), legend.text=element_text(size=7),panel.spacing=unit(0,"lines"), strip.background=element_rect(fill="white",color=NA), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(size=0.2), panel.spacing.y=unit(0,"lines"), panel.spacing.x=unit(1,"lines"))+xlab("Odds ratio of sex-biased expression compared to protein-coding genes")+guides(color=guide_legend(title="Gene type"))+scale_color_manual(values=c("#4daf4a","#ff7f00","#f781bf"))+scale_linetype_manual(values=c("solid","solid"))
dev.off()


