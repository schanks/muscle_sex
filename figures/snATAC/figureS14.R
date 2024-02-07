library(data.table)
library(ggplot2)


dat=fread("15state_nocount_female.tab")
dat$cell=gsub("_"," ",dat$cell)
dat$UB[which(dat$UB=="Inf")]=NA
dat$LB[which(dat$LB==0)]=NA
dat$State=unlist(strsplit(dat$State, "_"))[seq(2, nrow(dat)*2, by=2)]
dat$Adjustment=rep("No adjustment", nrow(dat))
d=fread("15state_female.tab")
d$cell=gsub("_"," ",d$cell)
d$UB[which(d$UB=="Inf")]=NA
d$LB[which(d$LB==0)]=NA
d$State=unlist(strsplit(d$State, "_"))[seq(2, nrow(d)*2, by=2)]
d$Adjustment=rep("Peak count adjusted", nrow(d))
dat=rbind(dat, d)
dat$Annotation=rep("E108: Female", nrow(dat))
f=dat
dat=fread("15state_nocount_male.tab")
dat$cell=gsub("_"," ",dat$cell)
dat$UB[which(dat$UB=="Inf")]=NA
dat$LB[which(dat$LB==0)]=NA
dat$State=unlist(strsplit(dat$State, "_"))[seq(2, nrow(dat)*2, by=2)]
dat$Adjustment=rep("No adjustment", nrow(dat))
d=fread("15state_male.tab")
d$cell=gsub("_"," ",d$cell)
d$UB[which(d$UB=="Inf")]=NA
d$LB[which(d$LB==0)]=NA
d$State=unlist(strsplit(d$State, "_"))[seq(2, nrow(d)*2, by=2)]
d$Adjustment=rep("Peak count adjusted", nrow(d))
dat=rbind(dat, d)
dat$Annotation=rep("E107: Male", nrow(dat))
m=dat
dat=fread("15state_nocount_consensus.tab")
dat$cell=gsub("_"," ",dat$cell)
dat$UB[which(dat$UB=="Inf")]=NA
dat$LB[which(dat$LB==0)]=NA
dat$State=unlist(strsplit(dat$State, "_"))[seq(2, nrow(dat)*2, by=2)]
dat$Adjustment=rep("No adjustment", nrow(dat))
d=fread("15state_consensus.tab")
d$cell=gsub("_"," ",d$cell)
d$UB[which(d$UB=="Inf")]=NA
d$LB[which(d$LB==0)]=NA
d$State=unlist(strsplit(d$State, "_"))[seq(2, nrow(d)*2, by=2)]
d$Adjustment=rep("Peak count adjusted", nrow(d))
dat=rbind(dat, d)
dat$Annotation=rep("Consensus", nrow(dat))
dat=rbind(dat, m, f)
dat$yval=rep(1, nrow(dat))
dat$yval[which(dat$Annotation=="Consensus")]=3
dat$yval[which(dat$Annotation=="E107: Male")]=2
dat$yval[which(dat$Adjustment=="No adjustment")]=dat$yval[which(dat$Adjustment=="No adjustment")]+0.3


tiff("~/plot.tiff", units="mm",height=180, width=180, res=300)
ggplot(dat, aes(x=OR, y=yval, color=Adjustment))+theme_bw()+geom_vline(xintercept=1, linetype="longdash")+geom_point(aes(shape=Annotation), size=2)+geom_segment(aes(yend=yval,x=LB, xend=UB), linewidth=1)+scale_x_log10(position="top")+scale_alpha_manual()+scale_shape_manual(values=c(19,8,25))+scale_color_manual(values=c("#1b9e77","#e7298a"))+facet_grid(State~cell, switch="both")+theme(axis.text.x=element_text(size=7), strip.text.y.left=element_text(size=7, angle=0),strip.text.x=element_text(size=7),axis.title.x=element_text(size=7),legend.title=element_text(size=7), legend.text=element_text(size=7),panel.spacing.x=unit(1,"lines"), panel.spacing.y=unit(0, "lines"),strip.background=element_rect(fill="white",color=NA), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(size=0.3,color="black"))+xlab("Odds ratio of sex-biased accessibility compared to Quiescent/Low Signal")
dev.off()


