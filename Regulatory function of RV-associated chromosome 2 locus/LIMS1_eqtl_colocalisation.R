rm(list = ls(all=TRUE))

library(coloc)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#read in data
#gwas_beta/gwas_se: association statistics for the Uganda anti-RV IgA GWAS
#mono_beta/mono_se: LIMS1 cis eQTL mapping statistics in monocytes from: https://pubmed.ncbi.nlm.nih.gov/24604202/
#bcell_beta/bcell_se: LIMS1 cis eQTL mapping statistics in B cells from: https://pubmed.ncbi.nlm.nih.gov/22446964/
#nk_beta/nk_se: LIMS1 cis eQTL mapping statistics in NK cells from: https://pubmed.ncbi.nlm.nih.gov/35835762/
#cd4_beta/cd4_se: LIMS1 cis eQTL mapping statistics in CD4+ T cells from: https://pubmed.ncbi.nlm.nih.gov/28248954/
#cd8_beta/cd8_se: LIMS1 cis eQTL mapping statistics in CD8+ T cells from: https://pubmed.ncbi.nlm.nih.gov/28248954/
#pbmc_beta/pbmc_se: LIMS1 cis eQTL mapping statistics in PBMCs from: https://pubmed.ncbi.nlm.nih.gov/39056362/
total.data.out <- read.table("coloc_eqtl_summ_stats.txt", header = T)

#format into list format for coloc
gwas.list <- list(beta = total.data.out$gwas_beta, varbeta = (total.data.out$gwas_se)^2, type = "quant", N = 616, sdY = 1, snp = as.character(total.data.out$snp))
mono.list <- list(beta = total.data.out$mono_beta[which(!is.na(total.data.out$mono_beta))], varbeta = (total.data.out$mono_se[which(!is.na(total.data.out$mono_beta))])^2, type = "quant", N =400, sdY = 1, snp = as.character(total.data.out$snp[which(!is.na(total.data.out$mono_beta))]))
nk.list <- list(beta = total.data.out$nk_beta[which(!is.na(total.data.out$nk_beta))], varbeta = (total.data.out$nk_se[which(!is.na(total.data.out$nk_beta))])^2, type = "quant", N = 245, sdY = 1, snp = as.character(total.data.out$snp[which(!is.na(total.data.out$nk_beta))]))
bcell.list <- list(beta = total.data.out$bcell_beta[which(!is.na(total.data.out$bcell_beta))], varbeta = (total.data.out$bcell_se[which(!is.na(total.data.out$bcell_beta))])^2, type = "quant", N = 281, sdY = 1, snp = as.character(total.data.out$snp[which(!is.na(total.data.out$bcell_beta))]))
cd4.list <- list(beta = total.data.out$cd4_beta[which(!is.na(total.data.out$cd4_beta))], varbeta = (total.data.out$cd4_se[which(!is.na(total.data.out$cd4_beta))])^2, type = "quant", N = 293, sdY = 1, snp = as.character(total.data.out$snp[which(!is.na(total.data.out$cd4_beta))]))
cd8.list <- list(beta = total.data.out$cd8_beta[which(!is.na(total.data.out$cd8_beta))], varbeta = (total.data.out$cd8_se[which(!is.na(total.data.out$cd8_beta))])^2, type = "quant", N = 283, sdY = 1, snp = as.character(total.data.out$snp[which(!is.na(total.data.out$cd8_beta))]))
pbmc.list <- list(beta = total.data.out$pbmc_beta[which(!is.na(total.data.out$pbmc_beta))], varbeta = (total.data.out$pbmc_se[which(!is.na(total.data.out$pbmc_beta))])^2, type = "quant", N = 1012, sdY = 1, snp = as.character(total.data.out$snp[which(!is.na(total.data.out$pbmc_beta))]))

#assess evidence for colocalisation between RV immunity GWAS and LIMS cis eQTLs across immune cells subsets
coloc.abf(gwas.list, mono.list)
coloc.abf(gwas.list, nk.list)
coloc.abf(gwas.list, bcell.list)
coloc.abf(gwas.list, cd4.list)
coloc.abf(gwas.list, cd8.list)
coloc.abf(gwas.list, pbmc.list)



