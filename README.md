# RV_GWAS

**Host genetic determinants of rotavirus disease and immunity in African children.**

Code and source data for GWAS of rotavirus disease and immunity in Kenya, Uganda and Zambia.

**Abstract**
Rotavirus (RV) is a major cause of diarrhoeal disease and mortality in African children. The impact of RV vaccination on child health has been limited by low vaccine effectiveness in resource-limited settings. We urgently need better vaccines to control RV disease, delivery of which will require an improved understanding of RV biology. We hypothesised that the genetics of RV disease and immunity are likely to be mutually informative of RV biology. We performed genome-wide association studies of infection-induced and vaccine-derived RV immunity, hospitalised RV gastroenteritis (RVGE) risk and all-cause diarrhoeal disease risk in >3,000 Ugandan, Kenyan and Zambian children. In doing so, we expand on the role of human blood group antigen secretion in RV infection and immunity, resolving complex interactions between secretor status, ABO blood group, Lewis antigen and viral genotype. We further identify a separate, novel locus at chromosome 2q12, rs59241810:T, associated with increased anti-RV IgA responses and lower risk of RVGE and all-cause diarrhoea with dehydration. Variation at rs59241810 regulates monocyte LIMS1 expression, thereby impacting IFNy-induced hepatocyte growth factor (HGF) secretion, suggesting a role for monocyte-derived HGF in the generation of anti-RV IgA responses. Our data reveal the complex interplay between immunity, disease, pathogen variation and host genetics in RV infections, and identify LIMS1/HGF signalling as a target to design more effective anti-RV vaccines.

**Overview of repository**
* Rotavirus disease & immunity genome-wide association analysis:
  1. Scripts used to perform association for GWAS of anti-RV immunity in Uganda and RV disease risk in Kenya (RV_GWAS.sh). Association analysis was performed using SNPTEST v2.5.6: [Marchini J, *et al*. "A new multipoint method for genome-wide association studies by imputation of genotypes." *Nature Genetics.* 2009.](https://doi.org/10.1038/ng2088)
  2. Scripts to plot manhattan and QQ plots for GWAS of anti-RV immunity in Uganda (uganda_GWAS_manhs_qq.R) and RV disease risk in Kenya (kenya_GWAS_manhs_qq.R) from thinned summary statistics (uganda_GWAS_summ_stats.thinned.txt.gz & kenya_GWAS_summ_stats.thinned.txt.gz).
  3. Script to plot regional association plot for secretor region (secretor_RA_plot.R) using Uganda and Kenya GWAS summary statistics (secretor_RA.ld.txt), local gene coordinates (secretor_genes.txt) and local recombination rate (genetic_map_chr19_b37.txt).
  4. Script to plot regional association plot for associated chromosome 2 region (chr2_locus_RA_plot.R) using Uganda GWAS summary statistics (chr2_RA.ld.txt), local gene coordinates (chr2_genes.txt) and local recombination rate (genetic_map_chr2_b37.txt).
 
* Secretor status, ABO & Lewis:
  1. Script (secretor_lewis_ABO.R) and data (kenya_abo_sec_lew.txt; uganda_sec_lew.txt) to reproduce analysis exploring the association between secretor status, Lewis antigen status, ABO blood group and rotavirus disease risk and immunity.

  * Effect of chromosome 2 locus on RV immunity & disease:
  1. Script (secretor_lewis_ABO.R) and data (kenya_abo_sec_lew.txt; uganda_sec_lew.txt) to reproduce analysis exploring the association between the novel chromosome 2 locus and rotavirus disease risk, all-cause diarrhoeal disease risk and immunity.
