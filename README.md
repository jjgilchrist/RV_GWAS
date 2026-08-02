# RV_GWAS

**Host genetic determinants of rotavirus disease and immunity in African children.**

Code and source data for GWAS of rotavirus disease and immunity in Kenya, Uganda and Zambia.

**Abstract**
Rotavirus (RV) is a major cause of diarrhoeal disease and mortality in African children. The impact of RV vaccination on child health has been limited by low vaccine effectiveness in resource-limited settings. We urgently need better vaccines to control RV disease, delivery of which will require an improved understanding of RV biology. We hypothesised that the genetics of RV disease and immunity are likely to be mutually informative of RV biology. We performed genome-wide association studies of infection-induced and vaccine-derived RV immunity, hospitalised RV gastroenteritis (RVGE) risk and all-cause diarrhoeal disease risk in >3,000 Ugandan, Kenyan and Zambian children. In doing so, we expand on the role of human blood group antigen secretion in RV infection and immunity, resolving complex interactions between secretor status, ABO blood group, Lewis antigen and viral genotype. We further identify a separate, novel locus at chromosome 2q12, rs59241810:T, associated with increased anti-RV IgA responses and lower risk of RVGE and all-cause diarrhoea with dehydration. Variation at rs59241810 regulates monocyte LIMS1 expression, thereby impacting IFNgamma-induced hepatocyte growth factor (HGF) secretion, suggesting a role for monocyte-derived HGF in the generation of anti-RV IgA responses. Our data reveal the complex interplay between immunity, disease, pathogen variation and host genetics in RV infections, and identify LIMS1/HGF signalling as a target to design more effective anti-RV vaccines.

**Overview of repository**
* Rotavirus disease & immunity genome-wide association analysis:
  1. Scripts to perform GWAS of anti-RV immunity in Uganda and RV disease risk in Kenya (RV_GWAS.sh). Association analysis was performed using SNPTEST v2.5.1: [Marchini J, *et al*. "A new multipoint method for genome-wide association studies by imputation of genotypes." *Nature Genetics.* 2009.](https://doi.org/10.1038/ng2088)
  2. Scripts to plot manhattan and QQ plots for GWAS of anti-RV immunity in Uganda (uganda_GWAS_manhs_qq.R) and RV disease risk in Kenya (kenya_GWAS_manhs_qq.R) from thinned summary statistics (uganda_GWAS_summ_stats.thinned.txt.gz & kenya_GWAS_summ_stats.thinned.txt.gz).
  3. Script to plot regional association plot for secretor region (secretor_RA_plot.R) using Uganda and Kenya GWAS summary statistics (secretor_RA.ld.txt), local gene coordinates (secretor_genes.txt) and local recombination rate (genetic_map_chr19_b37.txt).
  4. Script to plot regional association plot for associated chromosome 2 region (chr2_locus_RA_plot.R) using Uganda GWAS summary statistics (chr2_RA.ld.txt), local gene coordinates (chr2_genes.txt) and local recombination rate (genetic_map_chr2_b37.txt).
 
* Secretor status, ABO & Lewis:
  1. Script (secretor_lewis_ABO.R) and data (kenya_abo_sec_lew.txt; uganda_sec_lew.txt) to illustrate analysis exploring the association between secretor status, Lewis antigen status, ABO blood group and rotavirus disease risk and immunity.

* Effect of chromosome 2 locus on RV immunity & disease:
  1. Script (chr_2_locus_rep.R) and data (uganda_diarrhoea.txt; zambia_chr2_locus.txt; kenya_RV_disease_chr2.txt; uganda_IgA_chr2.txt) to illustrate analysis exploring the association between the novel chromosome 2 locus and rotavirus disease risk, all-cause diarrhoeal disease risk and immunity.
 
* Regulatory function of RV-associated chromosome 2 locus:
  1. Script (LIMS1_eQTL_colocalisation.R) and data (coloc_eqtl_summ_stats.txt) to reproduce colocalisation analysis effect of the RV-associated chromosome 2 locus on gene expression in immune cells. Colocalisation analysis was performed using Coloc v5.0.1: [Giambartolomei C, *et al*. "Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics." *PLoS Genetics.* 2014.](https://doi.org/10.1371/journal.pgen.1004383)
  2. Script (LIMS1_cytokine_correlation.R) and data (IFN_monocytes_cytokine.txt) to illustrate correlation between LIMS1 expression and secretion of HGF from stimulated (IFN-gamma 24 hours) monocytes.
 
**Data availability and sources**
Complete GWAS summary statistics will be made available through the NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics). Individual level genotype and phenotype data will be deposited with the European Genome-Phenome Archive.

Colocalisation analysis makes use of previously-published eQTL and cytokine QTL mapping data:
* T cells [Kasela S, *et al*. "Pathogenic implications for autoimmune mechanisms derived by comparative eQTL analysis of CD4+ versus CD8+ T cells." *PLoS Genetics.* 2017.](https://doi.org/10.1371/journal.pgen.1006643)
* Monocytes [Fairfax BP, *et al*. "Innate immune activity conditions the effect of regulatory variants upon monocyte gene expression." *Science*. 2014.](https://doi.org/10.1126/science.1246949)
* B cells [Fairfax BP, *et al*. "Genetics of gene expression in primary immune cells identifies cell type-specific master regulators and roles of HLA alleles." *Nature Genetics*. 2012.](https://doi.org/10.1038/ng.2205)
* NK cells [Gilchrist JJ, *et al*. "Natural Killer cells demonstrate distinct eQTL and transcriptome-wide disease associations, highlighting their role in autoimmunity." *Nature Communications.* 2022.](https://doi.org/10.1038/s41467-022-31626-4)
* PBMCs [Wen J, *et al*. "Gene Expression and Splicing QTL Analysis of Blood Cells in African American Participants from the Jackson Heart Study." *Genetics*. 2024.](https://doi.org/10.1093/genetics/iyae098)
* Cytokines [Gilchrist JJ, *et al*. "Genetic determinants of cytokine production in activated human monocytes." *MedRxiv*. 2026.](https://doi.org/10.64898/2026.05.08.26352736)
