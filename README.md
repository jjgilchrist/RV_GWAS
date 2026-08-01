# RV_GWAS

**Host genetic determinants of rotavirus disease and immunity in African children.**

Code and source data for GWAS of rotavirus disease and immunity in Kenya, Uganda and Zambia.

**Abstract**
Rotavirus (RV) is a major cause of diarrhoeal disease and mortality in African children. The impact of RV vaccination on child health has been limited by low vaccine effectiveness in resource-limited settings. We urgently need better vaccines to control RV disease, delivery of which will require an improved understanding of RV biology. We hypothesised that the genetics of RV disease and immunity are likely to be mutually informative of RV biology. We performed genome-wide association studies of infection-induced and vaccine-derived RV immunity, hospitalised RV gastroenteritis (RVGE) risk and all-cause diarrhoeal disease risk in >3,000 Ugandan, Kenyan and Zambian children. In doing so, we expand on the role of human blood group antigen secretion in RV infection and immunity, resolving complex interactions between secretor status, ABO blood group, Lewis antigen and viral genotype. We further identify a separate, novel locus at chromosome 2q12, rs59241810:T, associated with increased anti-RV IgA responses and lower risk of RVGE and all-cause diarrhoea with dehydration. Variation at rs59241810 regulates monocyte LIMS1 expression, thereby impacting IFNy-induced hepatocyte growth factor (HGF) secretion, suggesting a role for monocyte-derived HGF in the generation of anti-RV IgA responses. Our data reveal the complex interplay between immunity, disease, pathogen variation and host genetics in RV infections, and identify LIMS1/HGF signalling as a target to design more effective anti-RV vaccines.

**Overview of repository**
* Rotavirus disease & immunity genome-wide association analysis:
  1. Scripts used to perform association for GWAS of anti-RV immunity in Uganda and RV disease risk in Kenya (RV_GWAS.sh).
  2. Scripts to plot manhattan and QQ plots for GWAS of anti-RV immunity in Uganda (uganda_GWAS_manhs_qq.R) and RV disease risk in Kenya (kenya_GWAS_manhs_qq.R) from thinned summary statistics (uganda_GWAS_summ_stats.thinned.txt.gz & kenya_GWAS_summ_stats.thinned.txt.gz).
  Association analysis was performed using SNPTEST v2.5.6: [Marchini J, *et al*. "A new multipoint method for genome-wide association studies by imputation of genotypes." *Nature Genetics.* 2009.](https://doi.org/10.1038/ng2088) 
