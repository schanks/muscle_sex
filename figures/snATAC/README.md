# Sex differences in cell type level chromatin accessibility

##Secondary analyses
chromstate.R: Tested enrichment of differential accessibility by sex by chromatin state adjusting for bins of peak count using consensus states from the male and female reference samples

chromstate_female.R: Same as above using only female reference sample

chromstate_male.R: Same as above using only male reference sample

chromstate_nocount.R: Tested enrichment of differential accessibility by sex by chromatin state without adjusting for bins of peak count using consensus states from the male and female reference samples

chromstate_female_nocount.R: Same as above using only female reference sample

chromstate_male_nocount.R: Same as above using only male reference sample

tfbs_overlap.R: Annotate ATAC-seq peaks with overlapping TFBS

peak_to_tss.R: Annotate genes with number of differentially accessible peaks upstream of TSS

enrich_gene_peak.R: Calculate enrichment of differential expression by sex by number and direction of differentially accessible peaks upstream of TSS

## Figures
figure4.R: Figure 4. Sex differences in cell-type specific chromatin accessibility in human skeletal muscle

figureS13.R: Supplementary Figure 13. Comparison of differential accessibility by sex across muscle cell types

figureS14.R: Supplementary Figure 14. Association of chromatin state with differential accessibility by sex 

