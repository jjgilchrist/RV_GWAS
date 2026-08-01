#!/bin/sh

#genome-wide association analysis for rotavirus disease risk in Kenyan children - logistic regression, additive association model
#individual-level genotypes will be deposited with the European Genome-phenome Archive

for i in {1..22}
do
./snptest_v2.5.1 \
-data RV.chr"$i".maf1_r7.hwe10.gen.gz \
RV_geno.sample \
-o chr"$i".RV2.pc4.sex.snptest.gz \
-frequentist add \
-method newml \
-pheno cc \
-cov_names sex pc1 pc2 pc3 pc4 \
-exclude_samples ibd_pca.excl \
-log chr"$i".RV2.pc4.log
done

#genome-wide association analysis for anti-rotavirus IgA immunity in Ugandan children - linear regression, additive association model
#individual-level genotypes will be deposited with the European Genome-phenome Archive

for i in {1..22}
do
./snptest_v2.5.1 \
-data RV_uganda.chr"$i".maf1_r3.hwe10.gen.gz \
RV_uganda.sample \
-o chr"$i".RV.pc2.sex.age.snptest.gz \
-frequentist add \
-method score \
-pheno iga \
-quantile_normalise_phenotypes \
-cov_names age sex pc1 pc2 \
-exclude_samples RV_uganda.excl \
-log chr"$i".RV.pc2.sex.age.log
done


