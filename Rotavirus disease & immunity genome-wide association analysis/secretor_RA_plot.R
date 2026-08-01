rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(LDlinkR)
library(ggplotify)
library(patchwork)
library(ggbio)
library(data.table)
library(GenomicRanges)
library(biomaRt)

#read in data - coordinate info, Uganda (immunity) and Kenya (disease) GWAS p-values and r2 to lead SNP (rs601338)
#R2 from 1000G AFR super-population

total <- read.table("secretor_RA.ld.txt", header = T)

#put R2 values into bins
total$bin_r2 <- 1
total$bin_r2[which(total$r2>0.2 & total$r2 <= 0.4)] <- 2
total$bin_r2[which(total$r2>0.4 & total$r2 <= 0.6)] <- 3
total$bin_r2[which(total$r2>0.6 & total$r2 <= 0.8)] <- 4
total$bin_r2[which(total$r2>0.8)] <- 5

#highlight rs601338 (secretor locus)
total$annotate <- 0
total$annotate[c(which(total$rsid=="rs601338"))] <- 1

#set colout palette
cols3 <- brewer.pal(8,"Paired")
cols2 <- brewer.pal(11,"Spectral")
cols <- brewer.pal(8,"Set2")

#plot out associations for Ugandan GWAS
secretor_total <- ggplot(total, aes(x=pos, y=-log10(uganda.p))) + 
  xlim(49206674-100000, 49206674+100000) +
  geom_point(data=subset(total, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(total, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(total, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(total, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(total, bin_r2==5), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value") +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(total, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_x = 25000) +
  annotate("text", x = (49275000), y = 14, label = quote(r^2), hjust = 0.5, size = 2) +
  annotate("point", x = (49295000), y = 12.5, size = 2, colour = cols2[1]) +
  annotate("text", x = (49275000), y = 12.5, label = c("0.8-1.0"), hjust = 0.5, size = 2) +
  annotate("point", x = (49295000), y = 11, size = 2, colour = cols2[2]) +
  annotate("text", x = (49275000), y = 11, label = c("0.6-0.8"), hjust = 0.5, size = 2) +
  annotate("point", x = (49295000), y = 9.5, size = 2, colour = cols2[3]) +
  annotate("text", x = (49275000), y = 9.5, label = c("0.4-0.6"), hjust = 0.5, size = 2) +
  annotate("point", x = (49295000), y = 8, size = 2, colour = cols2[5]) +
  annotate("text", x = (49275000), y = 8, label = c("0.2-0.4"), hjust = 0.5, size = 2) +
  annotate("point", x = (49295000), y = 6.5, size = 2, colour = cols[8]) +
  annotate("text", x = (49275000), y = 6.5, label = c("<0.2"), hjust = 0.5, size = 2) +
annotate("text", x = (49107000), y = 14, label = c("AntiRV IgA"), hjust = 0, size=6)

#plot out associations for Kenyan GWAS
secretor_rep <- ggplot(total, aes(x=pos, y=-log10(kenya.p))) + 
  xlim(49206674-100000, 49206674+100000) +
  geom_point(data=subset(total, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(total, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(total, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(total, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(total, bin_r2==5), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value") +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(total, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_x = 25000) +
  annotate("text", x = (49107000), y = 7, label = c("RV disease"), hjust = 0, size=6)

#read in coordinates of genes in region (b37)
genes <- read.table("secretor_genes.txt", header = T)

#plot out gene positions
plot.range <- c(49206674-100000, 49206674+100000)
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(49206674-100000, 49206674+100000) +
  ylim(c(1.9,2.1)) +
  geom_segment(data = genes[seq(15,29,2),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(16,26,2),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(15,29,2),], aes(x=mid, label=external_gene_name), size=3, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(16,26,2),], aes(x=mid, label=external_gene_name), size=3, col = c("black"),
                   nudge_y =0.05, segment.color = NA) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#read in recombination rate (b37)
recomb <- read.table("genetic_map_chr19_b37.txt", header = F)
colnames(recomb) <- c("position", "rate", "map")
recomb$position <- as.numeric(recomb$position)
recomb$rate <- as.numeric(recomb$rate)
recomb_sec <- subset(recomb, position>49206674-100000 & position<49206674+100000)


recomb_rate <- ggplot(recomb_sec, aes(x=position, y=rate)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 19") +
  scale_x_continuous(breaks=c(49100000,49200000,49300000),
                     labels=c("49.1Mb", "49.2Mb", "49.3Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

#plot out association plots (2x), genes and recombinations rate
p1 <- (secretor_total/secretor_rep/genes.plot/recomb_rate) + plot_layout(heights = c(3, 3, 1, 1))


ggplot2::ggsave(
  "secretor_region_RA.jpg",
  width = 7,
  height = 5,
  dpi = 300
)



