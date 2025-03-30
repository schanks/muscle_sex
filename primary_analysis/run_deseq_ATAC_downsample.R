##!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

.libPaths("/net/snowwhite/home/aujackso/R/x86_64-pc-linux-gnu-library/4.0");
shhh(library(dplyr)) 
library(readr)
library(ggplot2)
shhh(library(DESeq2)) 
options(scipen=1, digits=5)

invnorm = function(x, seed) {
  set.seed(seed)
  qnorm((rank(x, na.last="keep", ties.method="random") - 0.5) / sum(!is.na(x)))
}


args=commandArgs(trailingOnly=TRUE)
## Usage: Rscript run_deseq.R i j
if (length(args)==0) {
  stop("Two arguments must be supplied (cell i trait j)\n", call.=FALSE)
} 
i<-as.numeric(args[1])
j<-as.numeric(args[2])

#  no Type_1 for downsampling  
celltype<-c("Adipocyte", "Endothelial", "Macrophage", "Mesenchymal_Stem_Cell", "Neuromuscular_junction", "Neuronal", "Satellite_Cell", "Smooth_Muscle", "T_cell", "Type_2a", "Type_2x")
pds.orig<-readRDS(paste0("~/sn_muscle_2023/data/ATAC.pds.orig.Type_1.Rds")) #280 samples
eds.orig<-read_tsv(paste0("/net/snowwhite/home/schanks/muscle/snatac/downsample/Type_1_",celltype[i],".downsample.tab"))

traits<-c("SEX")      
traits2<-c("SEX.M")   
	
eds.0<-read_tsv(paste0("/net/snowwhite/home/schanks/muscle/snatac/downsample/Type_1_",celltype[i],".downsample.tab"))
	peaks<-eds.0$peak
	eds <- eds.0 %>% dplyr::select(-peak) %>% as.matrix()
	rownames(eds)<-peaks
	print(dim(eds))

	who.use<-colnames(eds)

	pds.0 <- pds.orig #Type_1
	print(dim(pds.0))
	pds <- pds.0 %>% dplyr::slice(match(who.use, labelcode))
	print(dim(pds))
	
	if (traits[j] %in% c("SEX","AGE")) {
	 	pds<-pds %>% 
		dplyr::select(SEX, labelcode, AGE, area_name, tss_enrichment, batch, fraction_mitochondrial_median, tot_nuc_allcelltypes) %>% 
		filter(complete.cases(.)) 
	}
	
	if (!(traits[j] %in% c("SEX","AGE"))) {
	 	pds<-pds %>% 
		dplyr::select(traits[j], SEX, labelcode, AGE, area_name, tss_enrichment, batch, fraction_mitochondrial_median, tot_nuc_allcelltypes) %>% 
		filter(complete.cases(.)) 
	}
	
	pds<-pds %>% 
		mutate(invtssenrich=invnorm(tss_enrichment,1), invage=invnorm(AGE,1), 
		invmito=invnorm(fraction_mitochondrial_median,1), 
		invtotnuc=invnorm(tot_nuc_allcelltypes, 1))  %>%
		mutate(batch2=gsub("-NM","",batch)) %>% 
		mutate(batch3=paste0("B",as.character(batch2))) %>%
		dplyr::select(-batch, -batch2) %>% 
		dplyr::rename(batch=batch3)	



	if (traits[j]=="SEX") {
		fstr<-"~ as.factor(batch) + as.factor(area_name) + invage + invtssenrich + invmito + invtotnuc + as.factor(SEX)" 
		redstr<-"~ as.factor(batch) + as.factor(area_name) + invage + invtssenrich + invmito + invtotnuc" 
	}
	if (traits[j]=="T2D") {
		fstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invage + invtssenrich + invmito + invtotnuc + T2D"
		redstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invage + invtssenrich + invmito + invtotnuc"
	}
	if (traits[j]=="AGE") {
		pds$x<- pds %>% pull(traits[j])
		pds$x<-invnorm(pds$x, j)
		fstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invtssenrich + invmito + invtotnuc + x"
		redstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invtssenrich + invmito + invtotnuc"
	}
	if (!(traits[j] %in% c("SEX", "AGE", "T2D"))) {
		pds$x<- pds %>% pull(traits[j])
		pds$x<-invnorm(pds$x, j)
		fstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invage + invtssenrich + invmito + invtotnuc + x"
		redstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invage + invtssenrich + invmito + invtotnuc"
	}
	
	# beware count 0 samples
	totct<-colSums(eds) %>% as.data.frame()
	colnames(totct)<-"totct"
	totct<-tibble::rownames_to_column(totct,'labelcode')
	who.use<-totct %>% filter(totct>0) %>% pull(labelcode)
	eds.use<-eds %>% as.data.frame() %>% dplyr::select(all_of(who.use)) %>% as.matrix()
	pds.use<-pds %>% as.data.frame() %>% dplyr::slice(match(who.use, labelcode))
	dds <- DESeqDataSetFromMatrix(countData = eds.use,
            colData = pds.use, design = as.formula(fstr)) 
	dds<-scran::computeSumFactors(dds)   
	dds<-DESeq(dds, test="LRT",reduced=as.formula(redstr), useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, fitType='local')  
	
	print(resultsNames(dds))
	stats<-mcols(dds) %>% as.data.frame() 
	stats$peak=rownames(stats)
	stats$cell=celltype[i]
	stats$trait=traits[j]
	converge<-stats %>% dplyr::select(peak, fullBetaConv)
	
	res = as(results(dds), "data.frame") 
	res$peak=rownames(res)
	res$cell=celltype[i]
	res$trait=traits2[j]
	res$n=nrow(pds)
	res<-merge(res,converge,by="peak")

res1<-res %>% arrange(pvalue)

		
write.table(format(res1, digits=8), 
file=paste0("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.ATAC/downsampled/Type_1_",celltype[i],".",traits2[j],".results.tab"), sep="\t", quote=F, col.names=T, row.names=F)

