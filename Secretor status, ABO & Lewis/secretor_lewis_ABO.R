rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(nnet)
library(patchwork)

#set colours
cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(10,"Paired")

#read in secretor, lewis, ABO and phenotype data (Kenya) 
#pc1-4: genotype principal components
#cc: case-control status (1 = admission with RV diarrhoea)
#p_geno: RV p genotype
#rsid.19.49206674.G.A.rec: secretor status (1 = non-secretors)
#bld_grp: ABO blood group
#O_nonO: O vs nonO (A, B, AB) blood group
#lewis_neg: Lewis antigen status (1 = Lewis negative)
ken.d <- read.table("/Users/jamesgilchrist/Documents/rotavirus/RV_CIDRZ_project/github_rep/secretor_abo_lewis/kenya_abo_sec_lew.txt", header = T)

#secretor status is associated with overall risk of hospital admission with rotavirus diarrhoea
sec.ken.total <- glm(cc~rsid.19.49206674.G.A.rec+sex+pc1+pc2+pc3+pc4, data = ken.d, family = "binomial")
summary(sec.ken.total)

#multinomial logistic regression models to define RV P genotypes associated with secretor status

#first set control samples as reference level
ken.d$p_geno <- factor(ken.d$p_geno)
ken.d$p_geno <- relevel(ken.d$p_geno, ref = "control")

#secretor status is associated with P[4] and P[8] genotype disease but not P[6]
sec.ken.geno <- multinom(p_geno ~ rsid.19.49206674.G.A.rec + sex +pc1+pc2+pc3+pc4, data=ken.d)
zvalues <- summary(sec.ken.geno)$coefficients / summary(sec.ken.geno)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#plot out association between secretor status, RV genotype and disease risk
#first extract betas and standard errors from models
label <- c("Total", "P4", "P6", "P8")
mean  <- c(summary(sec.ken.total)$coef[2,1],
           summary(sec.ken.geno)$coefficients[1,2],
           summary(sec.ken.geno)$coefficients[2,2],
           summary(sec.ken.geno)$coefficients[3,2])
se  <- c(summary(sec.ken.total)$coef[2,2],
         summary(sec.ken.geno)$standard.errors[1,2],
         summary(sec.ken.geno)$standard.errors[2,2],
         summary(sec.ken.geno)$standard.errors[3,2])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- "rs601338:AA"

df <- data.frame(label, mean, lower, upper, facet)

#plot forest plot of effects of secretor status
sec.geno.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
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
sec.geno.fp

#next explore effect of ABO blood group
#set blood group O as baseline
ken.d$bld_grp <- factor(ken.d$bld_grp, levels=c("O", "A", "B", "AB"))

#ABO blood group is associated with risk of hospitalised diarrhoea in Kenya
lm1 <- glm(cc~sex+pc1+pc2+pc3+pc4+bld_grp, data = ken.d, family = "binomial")
summary(lm1)

#ABO blood group associattion is restricted to risk of P[8] genotype disease
m1 <- multinom(p_geno ~ bld_grp + sex +pc1+pc2+pc3+pc4, data=ken.d)
zvalues <- summary(m1)$coefficients / summary(m1)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#evidence of interaction between blood group and secretor status in P[8] genotype disease risk
m1 <- multinom(p_geno ~ rsid.19.49206674.G.A.rec*O_nonO + sex +pc1+pc2+pc3+pc4, data=ken.d)
zvalues <- summary(m1)$coefficients / summary(m1)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#plot out effect of blood group on RV disease risk stratified by secretor status
#estimate effect of blood group on overall disease risk among secretors
total.sec.ests <- glm(cc ~ bld_grp + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==0), family = "binomial")
summary(total.sec.ests)
#estimate effect of blood group on overall disease risk among non-secretors
total.nonsec.ests <- glm(cc ~ bld_grp + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==1), family = "binomial")
summary(total.nonsec.ests)
#estimate effect of blood group on RV genotype-specific disease risk among secretors
sec.geno.ests <- multinom(p_geno ~ bld_grp + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==0))
zvalues <- summary(sec.geno.ests)$coefficients / summary(sec.geno.ests)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2
#estimate effect of blood group on per RV genotype-specific disease risk among non-secretors
nonsec.geno.ests <- multinom(p_geno ~ bld_grp + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==1))
zvalues <- summary(nonsec.geno.ests)$coefficients / summary(nonsec.geno.ests)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#collect betas and standard error into matrix
label <- rep(c("Total", "P4", "P6", "P8"),6)
mean  <- c(summary(total.sec.ests)$coefficients[2,1],
           summary(sec.geno.ests)$coefficients[1,2],
           summary(sec.geno.ests)$coefficients[2,2],
           summary(sec.geno.ests)$coefficients[3,2],
           summary(total.nonsec.ests)$coefficients[2,1],
           summary(nonsec.geno.ests)$coefficients[1,2],
           summary(nonsec.geno.ests)$coefficients[2,2],
           summary(nonsec.geno.ests)$coefficients[3,2],
           summary(total.sec.ests)$coefficients[3,1],
           summary(sec.geno.ests)$coefficients[1,3],
           summary(sec.geno.ests)$coefficients[2,3],
           summary(sec.geno.ests)$coefficients[3,3],
           summary(total.nonsec.ests)$coefficients[3,1],
           NA,
           summary(nonsec.geno.ests)$coefficients[2,3],
           summary(nonsec.geno.ests)$coefficients[3,3],
           summary(total.sec.ests)$coefficients[4,1],
           summary(sec.geno.ests)$coefficients[1,4],
           summary(sec.geno.ests)$coefficients[2,4],
           summary(sec.geno.ests)$coefficients[3,4],
           NA,
           NA,
           NA,
           NA)

se  <- c(summary(total.sec.ests)$coefficients[2,1],
         summary(sec.geno.ests)$standard.errors[1,2],
         summary(sec.geno.ests)$standard.errors[2,2],
         summary(sec.geno.ests)$standard.errors[3,2],
         summary(total.nonsec.ests)$coefficients[2,1],
         summary(nonsec.geno.ests)$standard.errors[1,2],
         summary(nonsec.geno.ests)$standard.errors[2,2],
         summary(nonsec.geno.ests)$standard.errors[3,2],
         summary(total.sec.ests)$coefficients[3,1],
         summary(sec.geno.ests)$standard.errors[1,3],
         summary(sec.geno.ests)$standard.errors[2,3],
         summary(sec.geno.ests)$standard.errors[3,3],
         summary(total.nonsec.ests)$coefficients[3,1],
         NA,
         summary(nonsec.geno.ests)$standard.errors[2,3],
         summary(nonsec.geno.ests)$standard.errors[3,3],
         summary(total.sec.ests)$coefficients[4,1],
         summary(sec.geno.ests)$standard.errors[1,4],
         summary(sec.geno.ests)$standard.errors[2,4],
         summary(sec.geno.ests)$standard.errors[3,4],
         NA,
         NA,
         NA,
         NA)
lower <- mean-1.96*se
upper <- mean+1.96*se
facet1 <- c(rep(c(rep("Secretors",4), rep("Non-secretors",4)),3))
facet2 <- c(rep("A vs O",8), rep("B vs O",8), rep("AB vs O",8))
df <- data.frame(label, mean, lower, upper, facet1, facet2)

df$facet1 <- factor(df$facet1, levels = c("Secretors", "Non-secretors"))
df$facet2 <- factor(df$facet2, levels = c("A vs O", "B vs O", "AB vs O"))

#plot out forest plot
abo_sec_geno_int.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=1.5, size=2) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(8,8,2,8)])) + scale_colour_manual(values = c(cols[c(8,8,2,8)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + xlab("RV genotype") + scale_x_discrete(expand = expand_scale(add = 1)) +
  #scale_y_continuous(breaks=c(-0.1,0,0.1),
  #labels=c("-0.1", "0", "0.1")) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_grid(facet1~facet2) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text = element_text(
    size = 12, color = "white", face = "bold"))
abo_sec_geno_int.fp

#effect of lewis antigen status on overall RV disease risk in Kenyan children
total.lewis.ken <- glm(cc~lewis_neg+sex+pc1+pc2+pc3+pc4, data = ken.d, family = "binomial")
summary(total.lewis.ken)

#effect of lewis antigen status on RV genotype-specific disease riski n Kenyan children
geno.lewis.ken <- multinom(p_geno ~ lewis_neg + sex +pc1+pc2+pc3+pc4, data=ken.d)
zvalues <- summary(geno.lewis.ken)$coefficients / summary(geno.lewis.ken)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#compile effects and standard errors
label <- c("Total", "P4", "P6", "P8")
mean  <- c(summary(total.lewis.ken)$coef[2,1],
           summary(geno.lewis.ken)$coefficients[1,2],
           summary(geno.lewis.ken)$coefficients[2,2],
           summary(geno.lewis.ken)$coefficients[3,2])
se  <- c(summary(total.lewis.ken)$coef[2,2],
         summary(geno.lewis.ken)$standard.errors[1,2],
         summary(geno.lewis.ken)$standard.errors[2,2],
         summary(geno.lewis.ken)$standard.errors[3,2])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- c(rep("Le(-) vs Le(+)",4))

df <- data.frame(label, mean, lower, upper, facet)

df$facet <- factor(df$facet, levels = c("Le(-) vs Le(+)"))

#forest plot of effect of Lewis antigen status on RV disease risk overall and per RV genotype in Kenyan children
lewis.geno.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=1.5, size=2) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(4:1)])) + scale_colour_manual(values = c(cols[c(4:1)])) +
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
lewis.geno.fp


#Recapitulate this analysis stratifying by secretor status
#Effect of Lewis antigen status on total RV disease risk among secretors
total.sec.ests <- glm(cc ~ lewis_neg + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==0), family = "binomial")
summary(total.sec.ests)
#Effect of Lewis antigen status on total RV disease risk among non-secretors
total.nonsec.ests <- glm(cc ~ lewis_neg + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==1), family = "binomial")
summary(total.nonsec.ests)
#Effect of Lewis antigen status on genotype-specific RV disease risk among secretors
sec.geno.ests <- multinom(p_geno ~ lewis_neg + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==0))
zvalues <- summary(sec.geno.ests)$coefficients / summary(sec.geno.ests)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2
#Effect of Lewis antigen status on genotype-specific RV disease risk among non-secretors
nonsec.geno.ests <- multinom(p_geno ~ lewis_neg + sex +pc1+pc2+pc3+pc4, data=subset(ken.d, round(rsid.19.49206674.G.A.rec)==1))
zvalues <- summary(nonsec.geno.ests)$coefficients / summary(nonsec.geno.ests)$standard.errors
pnorm(abs(zvalues), lower.tail=FALSE)*2

#compile these estimates into a matrix
label <- rep(c("Total", "P4", "P6", "P8"),6)
mean  <- c(summary(total.sec.ests)$coefficients[2,1],
           summary(sec.geno.ests)$coefficients[1,2],
           summary(sec.geno.ests)$coefficients[2,2],
           summary(sec.geno.ests)$coefficients[3,2],
           summary(total.nonsec.ests)$coefficients[2,1],
           summary(nonsec.geno.ests)$coefficients[1,2],
           summary(nonsec.geno.ests)$coefficients[2,2],
           summary(nonsec.geno.ests)$coefficients[3,2])

se  <- c(summary(total.sec.ests)$coefficients[2,1],
         summary(sec.geno.ests)$standard.errors[1,2],
         summary(sec.geno.ests)$standard.errors[2,2],
         summary(sec.geno.ests)$standard.errors[3,2],
         summary(total.nonsec.ests)$coefficients[2,1],
         summary(nonsec.geno.ests)$standard.errors[1,2],
         summary(nonsec.geno.ests)$standard.errors[2,2],
         summary(nonsec.geno.ests)$standard.errors[3,2])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet1 <- c(rep("Secretors",4), rep("Non-secretors",4))
facet2 <- rep("Le(-) vs Le(+)",8)
df <- data.frame(label, mean, lower, upper, facet1, facet2)

df$facet1 <- factor(df$facet1, levels = c("Secretors", "Non-secretors"))
df$facet2 <- factor(df$facet2, levels = c("Le(-) vs Le(+)"))


#plot out effects of Lewis antigen status on RV genotype-specific disease risk stratified by secretor status
lewis_sec_geno_int.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=1.5, size=2) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(4:1)])) + scale_colour_manual(values = c(cols[c(4:1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + xlab("RV genotype") + scale_x_discrete(expand = expand_scale(add = 1)) +
  #scale_y_continuous(breaks=c(-0.1,0,0.1),
  #labels=c("-0.1", "0", "0.1")) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_grid(facet1~facet2) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text = element_text(
    size = 12, color = "white", face = "bold"))
lewis_sec_geno_int.fp

#explore the effect of lewis antigen status on anti-RV immunity in Ugandan children and its interaction with secretor status
#read in Lewis, secretor and phenotype data (Uganda)

#pc1-2: genotype principal components
#age: age at serum sampling (2 or 3 years)
#rsid.19.49206674.G.A.rec: secretor status (1 = non-secretors)
#lewis_neg: Lewis antigen status (1 = Lewis negative)
#norm.iga: inverse normal transformed serum anti-RV IgA levels
ug.d <- read.table("/Users/jamesgilchrist/Documents/rotavirus/RV_CIDRZ_project/github_rep/secretor_abo_lewis/uganda_sec_lew.txt", header = T)

#overall there is no effect of Lewis antigen status on anti-RV IgA responses in children
lm1 <- lm(norm.iga~lewis_neg+sex+age+pc1+pc2, data = ug.d)
summary(lm1)$coef

#plot this out as box plot
p_labels = data.frame(expt = c("base"), label = c("NS"))
box.lewis = ggplot(ug.d, aes(x=factor(round(lewis_neg,0)), y=norm.iga)) +
  geom_dotplot(binaxis="y", binwidth=0.04, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1), labels = c("Le(+)", "Le(-)")) + aes(fill = factor(round(lewis_neg,0)), col=factor(round(lewis_neg,0))) + scale_fill_manual(values = cols[c(1,1,1)]) + scale_colour_manual(values = cols[c(1,1,1)]) +
  ylab("Serum anti-RV IgA") + xlab("Lewis") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15), strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 15, color = "white", face = "bold")) + ylim(NA, 2.5) +
  geom_text(x=1.5, y=2.3, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)
box.lewis

#however there is evidence for a significant interaction between secretor status and lewis antigen
lm1 <- lm(norm.iga~lewis_neg*rsid.19.49206674.G.A.rec+sex+age+pc1+pc2, data = ug.d)
summary(lm1)$coef

#explore effect of Lewis antigen status on anti-RV IgA in secretors and non-secretors
lm1 <- lm(norm.iga~lewis_neg+age+pc1+pc2, data = subset(ug.d, round(rsid.19.49206674.G.A.rec,0)==1))
summary(lm1)

lm1 <- lm(norm.iga~lewis_neg+age+pc1+pc2, data = subset(ug.d, round(rsid.19.49206674.G.A.rec,0)==0))
summary(lm1)

#plot out these associations as box plots
ug.d$sec.status <- "Secretors"
ug.d$sec.status[which(round(ug.d$rsid.19.49206674.G.A.rec, 0)==1)] <- "Non-secretors"

ug.d$sec.status <- factor(ug.d$sec.status, levels = c("Secretors", "Non-secretors"))
p_labels = data.frame(sec.status = c("Secretors","Non-secretors"), label = c("italic(P)==0.07", "italic(P)==5%*%10^-9"))
p_labels$sec.status <- factor(p_labels$sec.status, levels = c("Secretors","Non-secretors"))

strat.lewis = ggplot(ug.d, aes(x=factor(round(lewis_neg,0)), y=norm.iga)) +
  geom_dotplot(binaxis="y", binwidth=0.04, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1), labels = c("Le(+)", "Le(-)")) + aes(fill = factor(sec.status), col=factor(sec.status)) + scale_fill_manual(values = cols2[c(1,7)]) + scale_colour_manual(values = cols2[c(1,7)]) +
  ylab("Serum anti-RV IgA") + xlab("Lewis") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15), strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 15, color = "white", face = "bold")) + ylim(NA, 2.5) +
  facet_wrap(~sec.status, ncol = 2) +
  geom_text(x=1.5, y=2.3, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)

strat.lewis

