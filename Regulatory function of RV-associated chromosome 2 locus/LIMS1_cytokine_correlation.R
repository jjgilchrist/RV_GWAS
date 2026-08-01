rm(list = ls(all=TRUE))

library(dplyr)
library(ggplot2)
library(RColorBrewer)

#set colour palette
cols <- brewer.pal(8,"Set2")

#read in data
#IFN_[cytokine]: [cytokine] secretion (normalised) from monocytes following 24 hours of IFN-gamma stimulation
#geno: rs59241810 genotype (additive)
#LIMS1.ifn: LIMS1 RNA expression (normalised) in monocytes following 24 hours of IFN-gamma stimulation
#RNA.pc1-pc25: first 25 principal components of gene expression data in monocytes following 24 hours of IFN-gamma stimulation
cytokine_RNA_LIMS1 <- read.table("~/Documents/rotavirus/RV_CIDRZ_project/github_rep/function/IFN_monocytes_cytokine.txt", header = T)

#calculate correlation between LIMS1 RNA expression and secretion of 28 cytokines from monocytes following 24 hours of IFN-gamma stimulation 
p.ifn<- c()
r.ifn <- c()
for (i in c(1:28)){
  p.ifn[i]<- summary(lm(cytokine_RNA_LIMS1[,i]~cytokine_RNA_LIMS1$LIMS1.ifn+factor(cytokine_RNA_LIMS1$sex)+cytokine_RNA_LIMS1$age+as.matrix(cytokine_RNA_LIMS1[,c(which(colnames(cytokine_RNA_LIMS1)=="RNA.pc1"):which(colnames(cytokine_RNA_LIMS1)=="RNA.pc25"))])))$coef[2,4]
  r.ifn[i] <- pcor.test(na.omit(cytokine_RNA_LIMS1)$LIMS1.ifn, na.omit(cytokine_RNA_LIMS1)[,i], as.matrix(na.omit(cytokine_RNA_LIMS1)[,c(29,30,33:57)]))$estimate
}

#collect results into matrix
out <- data.frame(cbind(colnames(cytokine_RNA_LIMS1)[c(1:28)], r.ifn, p.ifn))
#calculate FDR
out$fdr.ifn <- p.adjust(out$p.ifn, method = "fdr")
out$r.ifn <- as.numeric(out$r.ifn)
out$p.ifn <- as.numeric(out$p.ifn)
out$fdr.ifn <- as.numeric(out$fdr.ifn)

#print most signifacntly associated cytokine - only HGF secretion is significantly correlated wirt LIMS1 transcription
print(head(out[order(out$p.ifn),],10))

#plot out the correlation between HGF protein secretion and LIMS1 RNA expression
p_labels6.2 = data.frame(cyto = c("IFN - 24h"), 
                         label = c("italic(P)==2.7%*%10^-5"))
p_labels6a.2 = data.frame(cyto = c("IFN - 24h"), 
                          label = c("r==-0.24"))
lims1_cyto_ifn.hgf = ggplot(cytokine_RNA_LIMS1, aes(y=LIMS1.ifn, x=IFN_HGF)) +
  geom_point(col=cols[4], size = 2, alpha=0.7) +
  geom_smooth(method=lm, col=cols[8]) +
  ylab(expression(paste(italic("LIMS1"), " RNA expression"))) + xlab("HGF secretion") +
  facet_wrap(~cyto, ncol = 1) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=20),
                     axis.title=element_text(size=20), strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 20, color = "white", face = "bold")) + ylim(NA, 11.3) + #xlim(-2.5,NA) +
  geom_text(x=-2, y=11.2, aes(label=label), data=p_labels6.2, parse=TRUE, inherit.aes=F, size = 5, hjust = 0) +
  geom_text(x=-2, y=11, aes(label=label), data=p_labels6a.2, parse=TRUE, inherit.aes=F, size = 5, hjust = 0)
lims1_cyto_ifn.hgf

#then test is rs59241810 genotype (as a predictor of LIMS1 expression) is predictive of HGF secretion:

summary(lm(cytokine_RNA_LIMS1$IFN_HGF~cytokine_RNA_LIMS1$geno+
             factor(cytokine_RNA_LIMS1$sex)+
             cytokine_RNA_LIMS1$age))

#plot out forest plot of that association:
p_labels8 = data.frame(expt = c("base"), label = c("italic(P)==0.01"))
box.HGF.geno = ggplot(na.omit(cytokine_RNA_LIMS1), aes(x=factor(round(geno,0)), y=IFN_HGF)) +
  geom_dotplot(binaxis="y", binwidth=0.05, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1), labels = c("CC", "CT")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[c(4,4)]) + scale_colour_manual(values = cols[c(4,4)]) +
  ylab("HGF secretion") + xlab("rs59241810") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=20),
                     axis.title=element_text(size=20), strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 20, color = "white", face = "bold")) + #ylim(NA, 12) +
  geom_text(x=1.5, y=2, aes(label=label), data=p_labels8, parse=TRUE, inherit.aes=F, size = 8)
box.HGF.geno



