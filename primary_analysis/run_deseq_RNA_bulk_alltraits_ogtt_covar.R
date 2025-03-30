#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

.libPaths("/net/snowwhite/home/aujackso/R/x86_64-pc-linux-gnu-library/4.0");
shhh(library(dplyr)) 
library(readr)
library(ggplot2)
shhh(library(DESeq2)) 
library(edgeR)
library(stringr)
options(scipen=1, digits=5)

invnorm = function(x, seed) {
  set.seed(seed)
  qnorm((rank(x, na.last="keep", ties.method="random") - 0.5) / sum(!is.na(x)))
}


args=commandArgs(trailingOnly=TRUE)
## Usage: Rscript run_deseq.R j
if (length(args)==0) {
  stop("One argument must be supplied (trait j)\n", call.=FALSE)
} 
j<-as.numeric(args[1])


# h.traits doesn't include T2D or traits previously in bulk phenotype dataset from 6.1_bulk_pbulk_deseq_hg38.Rmd
h.traits<-c("S_Insu", "glu_fast_biopsy", "HEIGHT", "WEIGHT", "WAIST", "whr", "HIP", "S_ALAT", "S_hs_CRP", "fS_Kol", "fS_Kol_HDL", "fS_Kol_LDL_c", "fS_Trigly", "GL0", "GL30", "GL60", "GL120", "S_LipoA1", "S_LipoB", "fS_C_pept", "fS_Krea", "S_GT", "S_Uraat", "fS_C_pept_30", "S_Insu_30", "B_GHb_A1C", "glu_2h_biopsy", "daily_energy_expenditure", "total_daily_physical_activity", "moderate_physical_activity", "strenuous_physical_activity",  "p_insu", "p_insu_30", "p_insu_60", "p_insu_120", "RFM", "sbp", "dbp", "ApoB_A1_ratio", "Ins_AUC_0to30", "Glu_AUC_0to30", "InsSec30", "InsGenIn", "DIo", "CpepGenIn", "HOMA", "matsuda_4pt", "matsuda_3pt", "LTPA_duration_all", "LTPA_duration_cond", "LTPA_duration_noncond", "LTPA_energy_all", "LTPA_energy_cond", "LTPA_energy_noncond", "LTPA_energy_light", "LTPA_energy_moderate", "LTPA_energy_modvig", "LTPA_energy_vigorous", "sleep_biopsy", "totalcw", "sleep_24h","GL0")      


traits<-c("SEX", "age_at_biopsy", "T2D", "S_Insu", "bmi",  "glu_fast_biopsy", "HEIGHT", "WEIGHT", "WAIST", "whr", "HIP", "S_ALAT", "S_hs_CRP", "fS_Kol", "fS_Kol_HDL", "fS_Kol_LDL_c", "fS_Trigly", "GL0", "GL30", "GL60", "GL120", "S_LipoA1", "S_LipoB", "fS_C_pept", "fS_Krea", "S_GT", "S_Uraat", "fS_C_pept_30", "S_Insu_30", "B_GHb_A1C", "glu_2h_biopsy", "daily_energy_expenditure", "total_daily_physical_activity", "moderate_physical_activity", "strenuous_physical_activity",  "p_insu", "p_insu_30", "p_insu_60", "p_insu_120", "RFM", "sbp", "dbp", "ApoB_A1_ratio", "Ins_AUC_0to30", "Glu_AUC_0to30", "InsSec30", "InsGenIn", "DIo", "CpepGenIn", "HOMA", "matsuda_4pt", "matsuda_3pt", "LTPA_duration_all", "LTPA_duration_cond", "LTPA_duration_noncond", "LTPA_energy_all", "LTPA_energy_cond", "LTPA_energy_noncond", "LTPA_energy_light", "LTPA_energy_moderate", "LTPA_energy_modvig", "LTPA_energy_vigorous", "sleep_biopsy", "totalcw", "sleep_24h")      
traits2<-c("SEX.M", "AGE", "T2D", "INSFAST", "BMI", "GLUFAST", "HEIGHT", "WEIGHT", "WAIST", "WHR", "HIP", "ALAT", "CRP", "CHOL", "HDL", "LDL", "TG", "GL0", "GL30", "GL60", "GL120", "APOA1", "APOB", "CPEP", "CREAT", "GT", "URIC", "CPEP30", "S_INS30", "HBA1C", "GLU2H", "ENERGY", "PA_TOT", "PA_MOD", "PA_STREN", "INS0", "INS30", "INS60", "INS120", "RFM", "SBP", "DBP", "APOB.A1.RATIO", "INS_AUC_0_30", "GLU_AUC_0_30", "INSSEC30", "INSGENIN", "DI", "CPEPGENIN", "HOMA", "MATSUDA4", "MATSUDA3", "LTPA.DUR.ALL", "LTPA.DUR.COND", "LTPA.DUR.NONCOND", "LTPA.EN.ALL", "LTPA.EN.COND", "LTPA.EN.NONCOND", "LTPA.EN.LIGHT", "LTPA.EN.MOD", "LTPA.EN.MODVIG", "LTPA.EN.VIG", "SLEEP", "ALCOHOL", "SLEEP24H") 
 #65 traits 


p.bulk<-readRDS("~/sn_muscle_2023/data/phen.bulk.Rds") %>% 
		mutate(T2D=ifelse(ogtt=="T2D",1,ifelse(ogtt=="NGT",0,NA))) #115 T2D NAs
		
heather<-read_csv("~/sn_muscle_2023/data/tissue.csv") %>% select(labelcode, all_of(h.traits), ogtt_status_paper)

p.bulk<-merge(p.bulk, heather, by="labelcode")
  #lots of traits have missing values

# 22323 genes filtered for having at least 5 counts in at least 25% of the 268 samples in the bulk dataset after overlapping the samples with snRNA (pseudobulk)
#  Note that I am NOT dropping genes later based on 5 counts in at least 25% of samples with phenotypes!
bulk.filt<-readRDS("~/sn_muscle_2023/data/ct.bulk.filt.Rds") #22323 genes, 268 samples
annot.sn<-readRDS("~/sn_muscle_2023/data/annot.sn.Rds")



######################################################
## Base + cell type proportion covariates model - Filter for low counts     
# bulk counts ~  SEX + AGE + as.factor(mrna_batch) + RIN + TIN + as.factor(area_name) + GC + INSERT + ADI + ENDO + MACRO + MSC + NMJ + NEUR + SATEL + SM + TC + as.factor(ogtt_status_paper) + trait    
# Use Type 1 as reference cell type.     
######################################################

	eds<-bulk.filt

	if (traits[j]=="SEX") {
	pds <- p.bulk %>% 
		dplyr::select(sex, labelcode, age_at_biopsy, mrna_batch, area_name, rin, tin_median, GC_byRead_Mean, InsertSize_Median, Type_2a, Type_2x, Adipocyte, Endothelial, Macrophage, Mesenchymal_Stem_Cell, Neuromuscular_junction, Neuronal, Satellite_Cell, Smooth_Muscle, T_cell, ogtt_status_paper) 
		
		fstr<-"~  AGE + as.factor(mrna_batch) + RIN + TIN + as.factor(area_name) + GC + INSERT + ADI + ENDO + MACRO + MSC + NMJ + NEUR + SATEL + SM + TC + as.factor(ogtt_status_paper) + as.factor(sex)" 
	}
	if (traits[j]=="T2D") {
	pds <- p.bulk %>% 
		dplyr::select(T2D, sex, labelcode, age_at_biopsy, mrna_batch, area_name, rin, tin_median, GC_byRead_Mean, InsertSize_Median, Type_2a, Type_2x, Adipocyte, Endothelial, Macrophage, Mesenchymal_Stem_Cell, Neuromuscular_junction, Neuronal, Satellite_Cell, Smooth_Muscle, T_cell, ogtt_status_paper) %>%
		dplyr::filter(complete.cases(.)) 
		fstr<-"~  as.factor(sex) + AGE + as.factor(mrna_batch) + RIN + TIN + as.factor(area_name) + GC + INSERT + ADI + ENDO + MACRO + MSC + NMJ + NEUR + SATEL + SM + TC + as.factor(ogtt_status_paper) + T2D" 
	}
	if (traits[j]=="AGE") {
	pds <- p.bulk %>% 
		dplyr::select(sex, labelcode, age_at_biopsy, mrna_batch, area_name, rin, tin_median, GC_byRead_Mean, InsertSize_Median, Type_2a, Type_2x, Adipocyte, Endothelial, Macrophage, Mesenchymal_Stem_Cell, Neuromuscular_junction, Neuronal, Satellite_Cell, Smooth_Muscle, T_cell, ogtt_status_paper) 
		
		fstr<-"~  as.factor(sex) + as.factor(mrna_batch) + RIN + TIN + as.factor(area_name) + GC + INSERT + ADI + ENDO + MACRO + MSC + NMJ + NEUR + SATEL + SM + TC + as.factor(ogtt_status_paper) + AGE" 
	}
	if (!(traits[j] %in% c("SEX", "AGE", "T2D"))) {
	pds <- p.bulk %>% 
		dplyr::select(sex, labelcode, age_at_biopsy, mrna_batch, area_name, rin, tin_median, GC_byRead_Mean, InsertSize_Median, traits[j], Type_2a, Type_2x, Adipocyte, Endothelial, Macrophage, Mesenchymal_Stem_Cell, Neuromuscular_junction, Neuronal, Satellite_Cell, Smooth_Muscle, T_cell, ogtt_status_paper) %>%
		dplyr::filter(complete.cases(.)) 
		pds$x<- pds %>% pull(traits[j])
		pds$x<-invnorm(pds$x, j)
		fstr<-"~  as.factor(sex) + AGE + as.factor(mrna_batch) + RIN + TIN + as.factor(area_name) + GC + INSERT + ADI + ENDO + MACRO + MSC + NMJ + NEUR + SATEL + SM + TC + as.factor(ogtt_status_paper) + x" 
	}

#inverse normalize covariates based on available phenotype
pds <- pds %>% 
	mutate(AGE=invnorm(age_at_biopsy,1), RIN=invnorm(rin,1), TIN=invnorm(tin_median,1), GC=invnorm(GC_byRead_Mean,1), INSERT=invnorm(InsertSize_Median,1), TYPE2A = invnorm(Type_2a,3), TYPE2X=invnorm(Type_2x,3), ADI=invnorm(Adipocyte,3), ENDO=invnorm(Endothelial,3), MACRO=invnorm(Macrophage,3), MSC=invnorm(Mesenchymal_Stem_Cell,3), NMJ=invnorm(Neuromuscular_junction,3), NEUR=invnorm(Neuronal,3), SATEL=invnorm(Satellite_Cell,3), SM=invnorm(Smooth_Muscle,3), TC=invnorm(T_cell,3))

samps<-pds %>% pull(labelcode) #to select gene columns
eds<-as.data.frame(eds) %>% dplyr::select(all_of(samps)) %>% as.matrix()	

#  Using original DESeq bulk settings
	dds <- DESeqDataSetFromMatrix(countData = eds,
            colData = pds, design = as.formula(fstr)) 
	dds<-estimateSizeFactors(dds)
	dds<-estimateDispersions(dds) 
        dds<-nbinomWaldTest(dds, maxit=1000)
	stats<-mcols(dds) %>% as.data.frame() 
	stats$gene=rownames(stats)
	converge<-stats %>% dplyr::select(gene, betaConv)
	
	res = as(results(dds), "data.frame") 
	res$gene=rownames(res)
	res$data="bulk.basecelltype"
	res$trait=traits2[j]
	res$n=nrow(pds)
	res<-merge(res,converge,by="gene")
	
res1<-merge(res, annot.sn, by="gene") %>% arrange(pvalue)
write.table(format(res1, digits=8), 
file=paste0("~/sn_muscle_2023/output/DESeq.RNA/bulk_adjprop_ogtt_covar/hg38.filter22323.bulk.",traits2[j],".results.tab"), sep="\t", quote=F, col.names=T, row.names=F)




