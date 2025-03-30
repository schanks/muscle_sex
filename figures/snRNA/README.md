# Sex differences in cell type level gene expression

## Secondary analyses

t1_downsample.R: Downsample Type 1 muscle fiber gene counts and samples to approximately the same number in other cell types

genetype.R: Tested enrichment of differential expression by sex by gene type (protein-coding, lncRNA, pseudogene, other) adjusting for bins of gene count

genetype_aut.R: Same as above for autosomal genes only

genetype_nocount.R: Tested enrichment of differential expression by sex by gene type (protein-coding, lncRNA, pseudogene, other) without adjusting for gene count

genetype_nocount_aut.R: Same as above for autosomal genes only

rna_enrich_prep.R: Inverse normalization of signed -log10 p-values of differential expression by sex in muscle fiber types for RNA Enrich analysis


## Figures
figure2.R: Figure 2. Sex differences in cell type-specific gene expression in human skeletal muscle

figureS4.R: Supplementary Figure 4. Correlation between total read count (UMIs) and number of sex-biased genes detected across cell types
 
figureS5.R: Supplementary Figure 5. Comparison of differential expression by sex across muscle cell types 

figureS6.R: Supplementary Figure 6. Expression levels of LPP

figureS7.R: Supplementary Figure 7. Sensitivity analyses with and without adjusting for oral glucose tolerance test (OGTT) status

figureS8.R: Supplementary Figure 8. Association of gene type with sex-biased expression

figureS9.R: Supplementary Figure 9. Gene expression levels in muscle fibers of top 20 autosomal sex-biased genes in oxidative phosphorylation and caveola GO terms

