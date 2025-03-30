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
  
celltype<-c("Adipocyte", "Endothelial", "Macrophage", "Mesenchymal_Stem_Cell", "Neuromuscular_junction", "Neuronal", "Satellite_Cell", "Smooth_Muscle", "T_cell", "Type_1", "Type_2a", "Type_2x")
pds.orig<-readRDS(paste0("~/sn_muscle_2023/data/ATAC.pds.orig.",celltype[i],".Rds"))
eds.orig<-readRDS(paste0("~/sn_muscle_2023/data/ATAC.eds.orig.",celltype[i],".Rds"))
peaks<-read_tsv(paste0("~/sn_muscle_2023/data/atac_meanpkct_",celltype[i],".txt")) %>% filter(mean.pk.ct > 1) %>% pull(peak)
eds.orig<-eds.orig %>% filter(peak %in% all_of(peaks))

traits<-c("SEX", "AGE", "T2D", "INS", "BMI",  "GLU", "HEIGHT", "WEIGHT", "WAIST", "WHR", "HIP", "S_ALAT", "S_hs_CRP", "fS_Kol", "fS_Kol_HDL", "fS_Kol_LDL_c", "fS_Trigly", "GL0", "GL30", "GL60", "GL120", "S_LipoA1", "S_LipoB", "fS_C_pept", "fS_Krea", "S_GT", "S_Uraat", "fS_C_pept_30", "S_Insu_30", "B_GHb_A1C", "glu_2h_biopsy", "daily_energy_expenditure", "total_daily_physical_activity", "moderate_physical_activity", "strenuous_physical_activity",  "p_insu", "p_insu_30", "p_insu_60", "p_insu_120", "RFM", "sbp", "dbp", "ApoB_A1_ratio", "Ins_AUC_0to30", "Glu_AUC_0to30", "InsSec30", "InsGenIn", "DIo", "CpepGenIn", "HOMA", "matsuda_4pt", "matsuda_3pt", "LTPA_duration_all", "LTPA_duration_cond", "LTPA_duration_noncond", "LTPA_energy_all", "LTPA_energy_cond", "LTPA_energy_noncond", "LTPA_energy_light", "LTPA_energy_moderate", "LTPA_energy_modvig", "LTPA_energy_vigorous", "sleep_biopsy", "totalcw", "sleep_24h")      
traits2<-c("SEX.M", "AGE", "T2D", "INSFAST", "BMI", "GLUFAST", "HEIGHT", "WEIGHT", "WAIST", "WHR", "HIP", "ALAT", "CRP", "CHOL", "HDL", "LDL", "TG", "GL0", "GL30", "GL60", "GL120", "APOA1", "APOB", "CPEP", "CREAT", "GT", "URIC", "CPEP30", "S_INS30", "HBA1C", "GLU2H", "ENERGY", "PA_TOT", "PA_MOD", "PA_STREN", "INS0", "INS30", "INS60", "INS120", "RFM", "SBP", "DBP", "APOB.A1.RATIO", "INS_AUC_0_30", "GLU_AUC_0_30", "INSSEC30", "INSGENIN", "DI", "CPEPGENIN", "HOMA", "MATSUDA4", "MATSUDA3", "LTPA.DUR.ALL", "LTPA.DUR.COND", "LTPA.DUR.NONCOND", "LTPA.EN.ALL", "LTPA.EN.COND", "LTPA.EN.NONCOND", "LTPA.EN.LIGHT", "LTPA.EN.MOD", "LTPA.EN.MODVIG", "LTPA.EN.VIG", "SLEEP", "ALCOHOL", "SLEEP24H")   
	pds <- pds.orig %>% dplyr::filter(n_nuclei>=10) 
	# no samples should be dropped anymore, since dropping occurs in /analysis/3.1_deseq_snATAC.Rmd now
	
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
		mutate(batch3=paste0("B",as.character(batch))) %>%
		dplyr::select(-batch, -batch2) %>% 
		dplyr::rename(batch=batch3)	


	samps<-pds %>% pull(labelcode) #to select peak columns
	print(paste0("N samples ",celltype[i]," ",traits[j],": ", length(samps)) )
	
	eds<-eds.orig %>% dplyr::select(peak, all_of(samps))	
	peaks=eds$peak
	eds <- eds %>% dplyr::select(-peak) %>% as.matrix()
	rownames(eds)<-peaks
	
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
	dds <- DESeqDataSetFromMatrix(countData = eds,
            colData = pds, design = as.formula(fstr)) 
	dds<-scran::computeSumFactors(dds)   
	dds<-DESeq(dds, test="LRT",reduced=as.formula(redstr), useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, fitType="local")  
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
	
write.table(format(res, digits=8), 
	file=paste0("~/sn_muscle_2023/output/DESeq.ATAC/final_drop10nuc/results/",celltype[i],".",traits2[j],".results.tab"), sep="\t", quote=F, col.names=T, row.names=F)
	

