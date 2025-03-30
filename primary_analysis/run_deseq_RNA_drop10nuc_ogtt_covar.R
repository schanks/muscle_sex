# #!/usr/bin/env Rscript
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
  
pds.orig<-readRDS("~/sn_muscle_2023/data/pds.orig.Rds")
eds.orig<-readRDS("~/sn_muscle_2023/data/eds.orig.Rds")

#  Not analyzing Adipocyte, No T_cell (and no MFM as of 11/15/2022)  - but keep 12 non-MFM celltypes since that's how pds.orig and eds.orig were created
celltype<-c("Adipocyte", "Endothelial", "Macrophage", "Mesenchymal_Stem_Cell", "Neuromuscular_junction", "Neuronal", "Satellite_Cell", "Smooth_Muscle", "T_cell", "Type_1", "Type_2a", "Type_2x")

traits<-c("SEX", "AGE", "T2D", "INS", "BMI",  "GLU", "HEIGHT", "WEIGHT", "WAIST", "WHR", "HIP", "S_ALAT", "S_hs_CRP", "fS_Kol", "fS_Kol_HDL", "fS_Kol_LDL_c", "fS_Trigly", "GL0", "GL30", "GL60", "GL120", "S_LipoA1", "S_LipoB", "fS_C_pept", "fS_Krea", "S_GT", "S_Uraat", "fS_C_pept_30", "S_Insu_30", "B_GHb_A1C", "glu_2h_biopsy", "daily_energy_expenditure", "total_daily_physical_activity", "moderate_physical_activity", "strenuous_physical_activity",  "p_insu", "p_insu_30", "p_insu_60", "p_insu_120", "RFM", "sbp", "dbp", "ApoB_A1_ratio", "Ins_AUC_0to30", "Glu_AUC_0to30", "InsSec30", "InsGenIn", "DIo", "CpepGenIn", "HOMA", "matsuda_4pt", "matsuda_3pt", "LTPA_duration_all", "LTPA_duration_cond", "LTPA_duration_noncond", "LTPA_energy_all", "LTPA_energy_cond", "LTPA_energy_noncond", "LTPA_energy_light", "LTPA_energy_moderate", "LTPA_energy_modvig", "LTPA_energy_vigorous", "sleep_biopsy", "totalcw", "sleep_24h")      
traits2<-c("SEX.M", "AGE", "T2D", "INSFAST", "BMI", "GLUFAST", "HEIGHT", "WEIGHT", "WAIST", "WHR", "HIP", "ALAT", "CRP", "CHOL", "HDL", "LDL", "TG", "GL0", "GL30", "GL60", "GL120", "APOA1", "APOB", "CPEP", "CREAT", "GT", "URIC", "CPEP30", "S_INS30", "HBA1C", "GLU2H", "ENERGY", "PA_TOT", "PA_MOD", "PA_STREN", "INS0", "INS30", "INS60", "INS120", "RFM", "SBP", "DBP", "APOB.A1.RATIO", "INS_AUC_0_30", "GLU_AUC_0_30", "INSSEC30", "INSGENIN", "DI", "CPEPGENIN", "HOMA", "MATSUDA4", "MATSUDA3", "LTPA.DUR.ALL", "LTPA.DUR.COND", "LTPA.DUR.NONCOND", "LTPA.EN.ALL", "LTPA.EN.COND", "LTPA.EN.NONCOND", "LTPA.EN.LIGHT", "LTPA.EN.MOD", "LTPA.EN.MODVIG", "LTPA.EN.VIG", "SLEEP", "ALCOHOL", "SLEEP24H")   

	# Drop samples with <10 nuclei for the cell type
	pds <- pds.orig[[i]] %>%
		dplyr::filter(tot_nuc >= 10)  
	samps.drop10<-pds %>% pull(labelcode) #to select gene columns
	
	eds<-eds.orig[[i]] %>% dplyr::select(gene, all_of(samps.drop10))	

	#  Drop genes missing data on more than 75% of samples with >10 nuclei - keep genes with non0 counts for at least 25% of samples  
	gene<-eds$gene
	count.tmp<-eds %>% dplyr::select(-gene)
	nsamp<-ncol(count.tmp)
	count.tmp<-as_tibble(count.tmp) 
	ct.non0=rowSums(count.tmp!=0)
	prop.non0=ct.non0/nsamp
	keep.gene <- cbind.data.frame(gene, prop.non0) %>% dplyr::filter(prop.non0>=0.25) %>% pull(gene)
	eds<-eds %>% dplyr::slice(match(keep.gene, gene))
	print(paste0("N genes after excluding >=75% missing ",celltype[i],": ",nrow(eds)))

		
	#  Now drop samples based on available phenotype	
	print(dim(pds))
	if (traits[j] %in% c("SEX","AGE")) {
	 	pds<-pds %>% 
		dplyr::select(SEX, labelcode, AGE, med.frac.mito, area_name, batch, tot_nuc_allcelltypes, ogtt_status_paper) %>% 
		dplyr::filter(complete.cases(.)) 
	}
	if (!(traits[j] %in% c("SEX","AGE"))) {
	 	pds<-pds %>% 
		dplyr::select(traits[j], SEX, labelcode, AGE, med.frac.mito, area_name, batch, tot_nuc_allcelltypes, ogtt_status_paper)   %>% 
		dplyr::filter(complete.cases(.)) 
	}
	pds<-pds %>% mutate(invmito=invnorm(med.frac.mito,1), invage=invnorm(AGE,1), invnucall=invnorm(tot_nuc_allcelltypes,1)) 	

	samps<-pds %>% pull(labelcode) #to select gene columns
	print(paste0("N samples ",celltype[i]," ",traits[j],": ", length(samps)) )
	eds<-eds %>% dplyr::select(gene, all_of(samps))	
	genes=eds$gene
	eds <- eds %>% dplyr::select(-gene) %>% as.matrix()
	rownames(eds)<-genes

	if (traits[j]=="SEX") {
		fstr<-"~ as.factor(batch) + as.factor(area_name) + invage + invmito + invnucall + as.factor(ogtt_status_paper) + as.factor(SEX)" 
		redstr<-"~ as.factor(batch) + as.factor(area_name) + invage + invmito + invnucall + as.factor(ogtt_status_paper)" 
	}
	if (traits[j]=="T2D") {
		fstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX)+ invage + invmito + invnucall + as.factor(ogtt_status_paper) + T2D"
		redstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX)+ invage + invmito + invnucall + as.factor(ogtt_status_paper)"
	}
	if (traits[j]=="AGE") {
		pds$x<- pds %>% pull(traits[j])
		pds$x<-invnorm(pds$x, j)
		fstr<-"~ as.factor(batch) + as.factor(area_name) +  as.factor(SEX) + invmito + invnucall + as.factor(ogtt_status_paper) + x"
		redstr<-"~ as.factor(batch) + as.factor(area_name) +  as.factor(SEX) + invmito + invnucall + as.factor(ogtt_status_paper)"
	}
	if (!(traits[j] %in% c("SEX", "AGE", "T2D"))) {
		pds$x<- pds %>% pull(traits[j])
		pds$x<-invnorm(pds$x, j)
		fstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invage + invmito + invnucall + as.factor(ogtt_status_paper) + x"
		redstr<-"~ as.factor(batch) + as.factor(area_name) + as.factor(SEX) + invage + invmito + invnucall + as.factor(ogtt_status_paper)"
	}
	dds <- DESeqDataSetFromMatrix(countData = eds,
            colData = pds, design = as.formula(fstr)) 
	dds<-scran::computeSumFactors(dds)   
	dds<-DESeq(dds, test="LRT",reduced=as.formula(redstr), useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, fitType="local")  
	stats<-mcols(dds) %>% as.data.frame() 
	stats$gene=rownames(stats)
	stats$cell=celltype[i]
	stats$trait=traits[j]
	converge<-stats %>% dplyr::select(gene, fullBetaConv)
	
	res = as(results(dds), "data.frame") 
	res$gene=rownames(res)
	res$cell=celltype[i]
	res$trait=traits2[j]
	res$n=nrow(pds)
	res<-merge(res,converge,by="gene")
	
	
map<-read_tsv("~/sn_muscle_2023/data/gencode.v30.gene_lengths_tss.tsv") %>% select(gene_id, gene_name, chrom, gene_start, gene_end, gene_strand, gene_type)

res1<-merge(res, map, by.x="gene", by.y="gene_id", all.x=T) %>% arrange(pvalue)

		
write.table(format(res1, digits=8), 
file=paste0("~/sn_muscle_2023/output/DESeq.RNA/final_drop10nuc/results_ogtt_covar/",celltype[i],".",traits2[j],".results.tab"), sep="\t", quote=F, col.names=T, row.names=F)

