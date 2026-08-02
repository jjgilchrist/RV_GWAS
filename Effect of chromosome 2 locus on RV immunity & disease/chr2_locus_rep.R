rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(nnet)
library(patchwork)

#set colours
cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(10,"Paired")

#n.b. individual levels datasets used here are synthetic. They are provided to illustrate the work-flow and will not reproduce the study's results. Individual level genotype and phenotype data will be made available through the European Genome-phenome Archive.

#read in novel locus genotype (rs59241810) and phenotype (IgA) data (Uganda)
#pc1-2: genotype principal components
#age: age at serum sampling (2 or 3 years)
#rsid.2.109558305.C.T.add: rs59241810 genotype (additive)
#norm.iga: inverse normal transformed serum anti-RV IgA levels
ug.d <- read.table("uganda_IgA_chr2_synth.txt", header = T)

#Genotype at rs59241810 is associated with anti-RV IgA at 2/3 years in Ugandan children
summary(lm(norm.iga~rsid.2.109558305.C.T.add+pc1+pc2+sex+age, data = ug.d))

#plot our association (box plot)
p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.77%*%10^-12"))
box.novel = ggplot(ug.d, aes(x=factor(round(rsid.2.109558305.C.T.add,0)), y=norm.iga)) +
  geom_dotplot(binaxis="y", binwidth=0.025, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = factor(round(rsid.2.109558305.C.T.add,0)), col=factor(round(rsid.2.109558305.C.T.add,0))) + scale_fill_manual(values = cols[c(1,1,1)]) + scale_colour_manual(values = cols[c(1,1,1)]) +
  ylab("Serum anti-RV IgA") + xlab("rs59241810") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15), strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 15, color = "white", face = "bold")) + ylim(NA, 2.5) +
  geom_text(x=2, y=2.3, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)

box.novel

#test for effect of secretor status and rs59241810 genotype on anti-RV IgA responses in vaccinated Zambian children at 1 year of age
#read in data
#chr_2_locus: rs59241810 genotype (additive)
#secretor: secretor status (1 = non-secretors)
#feeding: infant feeding
#intervention: randomisation to additional vaccine dose
#norm.IgA: normalised serum anti-RV IgA at 12 months of age
zam.d <- read.table("zambia_chr2_locus_synth.txt", header = T)

#test for effect of rs59241810 and secretor status on anti-RV IgA levels in Zambia infants
summary(lm(norm.IgA~chr_2_locus+secretor+feeding+sex+intervention, data = zam.d))

#plot effect as box plot
zam.d$facet1 <- "Zambia"
p_labels = data.frame(expt = c("base"), label = c("italic(P)==0.049"))
box.chr2.zambia = ggplot(zam.d, aes(x=factor(chr_2_locus), y=norm.IgA)) +
  geom_dotplot(binaxis="y", binwidth=0.04, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = factor(chr_2_locus), col=factor(chr_2_locus)) + scale_fill_manual(values = cols[c(1,1,1)]) + scale_colour_manual(values = cols[c(1,1,1)]) +
  ylab("Serum anti-RV IgA") + xlab("rs59241810") +
  facet_wrap(~facet1, ncol = 1) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=20),
                     axis.title=element_text(size=25), strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 25, color = "white", face = "bold")) + ylim(NA, 2.7) +
  geom_text(x=2, y=2.5, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 8)
box.chr2.zambia

#test for an effect on disease risk at this locus in Kenyan children
#read in Kenyan data
#pc1-4: genotype principal components
#cc: case-control status (1 = admission with RV diarrhoea)
#p_geno: RV p genotype
#rsid.2.109558305.C.T.add: rs59241810 genotype (additive)
#age_cat: categorical variable - cases are grouped into "over" (age over 12 months at hospitalisation) and "under" (age under 12 months at hospitalisation)

ken.d <- read.table("kenya_RV_disease_chr2_synth.txt", header = T)

#no association between rs59241810 genotype and overall RV disease risk in Kenyan children
novel.ken.total <- glm(cc~rsid.2.109558305.C.T.add+sex+pc1+pc2+pc3+pc4, data = ken.d, family = "binomial")
summary(novel.ken.total)

#test for genotype-specific effect
#set P genotypes as factors with control as reference
ken.d$p_geno <- factor(ken.d$p_geno, levels = c("control", "[4]", "[6]", "[8]"))
ken.d$p_geno <- relevel(ken.d$p_geno, ref = "control")

#evidence of effect on P[8] genotype disease risk
novel.ken.geno <- multinom(p_geno ~ rsid.2.109558305.C.T.add + sex +pc1+pc2+pc3+pc4, data=ken.d)
zvalues <- summary(novel.ken.geno)$coefficients / summary(novel.ken.geno)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#compile effects into matrix
label <- c("Total", "P4", "P6", "P8")
mean  <- c(summary(novel.ken.total)$coef[2,1],
           summary(novel.ken.geno)$coefficients[1,2],
           summary(novel.ken.geno)$coefficients[2,2],
           summary(novel.ken.geno)$coefficients[3,2])
se  <- c(summary(novel.ken.total)$coef[2,2],
         summary(novel.ken.geno)$standard.errors[1,2],
         summary(novel.ken.geno)$standard.errors[2,2],
         summary(novel.ken.geno)$standard.errors[3,2])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- "rs59241810:T"
df <- data.frame(label, mean, lower, upper, facet)

#plot forest plot of effects of rs59241810 genotype on overall and per-genotype RV disease risk
novel.geno.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2, size=2) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(4:1)])) + scale_colour_manual(values = c(cols[c(4:1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + xlab("RV genotype") + scale_x_discrete(expand = expand_scale(add = 1)) +
  #scale_y_continuous(breaks=c(-0.1,0,0.1),
  #labels=c("-0.1", "0", "0.1")) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap(~facet, ncol = 1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
novel.geno.fp

#given that this effect appears to protect against disease through enhanced anti-RV immune reponses, would hypothesise that this effect would accrue with age
#significant effect of rs59241810 genotype carriage on overall disease risk seen in children over 12 months of age:
novel.ken.cat <- multinom(age_cat ~ rsid.2.109558305.C.T.add + sex +pc1+pc2+pc3+pc4, data=ken.d)
zvalues <- summary(novel.ken.cat)$coefficients / summary(novel.ken.cat)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#again this is driven by an effect on P[8] genotype disease
novel.ken.cat.young <- multinom(p_geno ~ rsid.2.109558305.C.T.add + sex +pc1+pc2+pc3+pc4, data=subset(ken.d,age_cat=="control" | age_cat=="under"))
zvalues <- summary(novel.ken.cat.young)$coefficients / summary(novel.ken.cat.young)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

novel.ken.cat.old <- multinom(p_geno ~ rsid.2.109558305.C.T.add + sex +pc1+pc2+pc3+pc4, data=subset(ken.d,age_cat=="control" | age_cat=="over"))
zvalues <- summary(novel.ken.cat.old)$coefficients / summary(novel.ken.cat.old)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#compile these effect into matrix
label <- c("Total", "P4", "P6", "P8")
mean  <- c(summary(novel.ken.cat)$coefficients[2,2],
           summary(novel.ken.cat.young)$coefficients[1,2],
           summary(novel.ken.cat.young)$coefficients[2,2],
           summary(novel.ken.cat.young)$coefficients[3,2],
           summary(novel.ken.cat)$coefficients[1,2],
           summary(novel.ken.cat.old)$coefficients[1,2],
           summary(novel.ken.cat.old)$coefficients[2,2],
           summary(novel.ken.cat.old)$coefficients[3,2]
           )
se  <- c(summary(novel.ken.cat)$standard.errors[2,2],
         summary(novel.ken.cat.young)$standard.errors[1,2],
         summary(novel.ken.cat.young)$standard.errors[2,2],
         summary(novel.ken.cat.young)$standard.errors[3,2],
         summary(novel.ken.cat)$standard.errors[1,2],
         summary(novel.ken.cat.old)$standard.errors[1,2],
         summary(novel.ken.cat.old)$standard.errors[2,2],
         summary(novel.ken.cat.old)$standard.errors[3,2])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- c(rep("<12 months",4),rep(">12 months",4))
df <- data.frame(label, mean, lower, upper, facet)


#plot out effects as forest plot
novel.geno.age.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2, size=2) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(4:1)])) + scale_colour_manual(values = c(cols[c(4:1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + xlab("RV genotype") + scale_x_discrete(expand = expand_scale(add = 1)) +
  #scale_y_continuous(breaks=c(-0.1,0,0.1),
  #labels=c("-0.1", "0", "0.1")) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap(~facet, ncol = 1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
novel.geno.age.fp

#test for an effect on all-cause diarrheoal disease risk at this locus in Ugandan children
#read in data - genotype and phenotype (diarrhoeal disease frequency) data for Ugandan children
#pc1-2: genotype principal components
#rsid.2.109558305.C.T.add: rs59241810 genotype (additive)
#totaldiaevent: total number reported episodes of diarrhoeal disease to age 10 years
#bin.dia.dehydr: ever had diarrhoea with dehydration
#bin.dia.dehydr: ever had diarrhoea with dysentery
#bin.dia.dehydr: ever had prolonged diarrhoea
ug2.d <- read.table("uganda_diarrhoea_synth.txt", header = T)


#is genotype at rs59241810 associated with all-cause diarrhoea?
#no evidence of association
lm1 <- glm(totaldiarevent~rsid.2.109558305.C.T.add+sex+pc1+pc2, data = ug2.d, family = poisson)
summary(lm1)

#is associated, however, with diarrhoea with dehydration
lm.dehyd <- glm(bin.dia.dehyr~rsid.2.109558305.C.T.add+sex+pc1+pc2, data = ug2.d, family = binomial)
summary(lm.dehyd)
#not with prolonged diarrhoea
lm.pro <- glm(bin.dia.pro~rsid.2.109558305.C.T.add+sex+pc1+pc2, data = ug2.d, family = binomial)
summary(lm.pro)
#not with dysentery
lm.dys <- glm(bin.dia.dys~rsid.2.109558305.C.T.add+sex+pc1+pc2, data = ug2.d, family = binomial)
summary(lm.dys)

#compile effects into matrix
label <- c("prolonged", "dehydration", "dysentery")
mean  <- c(summary(lm.pro)$coef[2,1],
           summary(lm.dehyd)$coef[2,1],
           summary(lm.dys)$coef[2,1])
se  <- c(summary(lm.pro)$coef[2,2],
         summary(lm.dehyd)$coef[2,2],
         summary(lm.dys)$coef[2,2])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- "rs59241810:T"
df <- data.frame(label, mean, lower, upper, facet)

#plot out effects of rs59241810 genotype on risk of diarrhoeal disease subtype

novel.geno.d.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2, size=2) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(1,8,8)])) + scale_colour_manual(values = c(cols[c(1,8,8)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + xlab("Diarrhoea subtype") + scale_x_discrete(expand = expand_scale(add = 1)) +
  #scale_y_continuous(breaks=c(-0.1,0,0.1),
  #labels=c("-0.1", "0", "0.1")) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=15),
                     axis.title.x=element_text(size=15), axis.title.y=element_blank()) +
  facet_wrap(~facet, ncol = 1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
novel.geno.d.fp



