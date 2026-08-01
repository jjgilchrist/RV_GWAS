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


#read in data - coordinate info, Uganda (immunity) GWAS p-values and r2 to lead SNP (rs59241810)
#R2 from 1000G AFR super-population

total <- read.table("chr2_RA.ld.txt", header = T)
total <- na.omit(total)

#put R2 values into bins
total$bin_r2 <- 1
total$bin_r2[which(total$r2>0.2 & total$r2 <= 0.4)] <- 2
total$bin_r2[which(total$r2>0.4 & total$r2 <= 0.6)] <- 3
total$bin_r2[which(total$r2>0.6 & total$r2 <= 0.8)] <- 4
total$bin_r2[which(total$r2>0.8)] <- 5

#highlight peak SNP (rs59241810)
total$annotate <- 0
total$annotate[c(which(total$rsid=="rs59241810"))] <- 1

#set colour palette
cols3 <- brewer.pal(8,"Paired")
cols2 <- brewer.pal(11,"Spectral")
cols <- brewer.pal(8,"Set2")

#plot out associations for Ugandan GWAS
chr2_total <- ggplot(total, aes(x=pos, y=-log10(uganda.p))) + 
  xlim(109558305-570000, 109558305+570000) +
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
  geom_text_repel( data=subset(total, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_x = 150000) +
  annotate("text", x = (109950000), y = 10, label = quote(r^2), hjust = 0.5, size = 2) +
  annotate("point", x = (110028305), y = 9.3, size = 2, colour = cols2[1]) +
  annotate("text", x = (109950000), y = 9.3, label = c("0.8-1.0"), hjust = 0.5, size = 2) +
  annotate("point", x = (110028305), y = 8.6, size = 2, colour = cols2[2]) +
  annotate("text", x = (109950000), y = 8.6, label = c("0.6-0.8"), hjust = 0.5, size = 2) +
  annotate("point", x = (110028305), y = 7.9, size = 2, colour = cols2[3]) +
  annotate("text", x = (109950000), y = 7.9, label = c("0.4-0.6"), hjust = 0.5, size = 2) +
  annotate("point", x = (110028305), y = 7.2, size = 2, colour = cols2[5]) +
  annotate("text", x = (109950000), y = 7.2, label = c("0.2-0.4"), hjust = 0.5, size = 2) +
  annotate("point", x = (110028305), y = 6.5, size = 2, colour = cols[8]) +
  annotate("text", x = (109950000), y = 6.5, label = c("<0.2"), hjust = 0.5, size = 2) +
annotate("text", x = (109078305), y = 11.5, label = c("AntiRV IgA"), hjust = 0, size=6)

#read in coordinates of genes in region (b37)
genes <- read.table("chr2_genes.txt", header = T)

#plot out gene positions
plot.range <- c(109558305-570000, 109558305+570000)
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(109558305-570000, 109558305+570000) +
  ylim(c(1.9,2.1)) +
  geom_segment(data = genes[seq(4,8,2),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(5,9,2),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(4,8,2),], aes(x=mid, label=external_gene_name), size=3, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(5,9,2),], aes(x=mid, label=external_gene_name), size=3, col = c("black"),
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
recomb <- read.table("genetic_map_chr2_b37.txt", header = F)
colnames(recomb) <- c("position", "rate", "map")
recomb$position <- as.numeric(recomb$position)
recomb$rate <- as.numeric(recomb$rate)
recomb_sec <- subset(recomb, position>109558305-570000 & position<109558305+570000)

recomb_rate <- ggplot(recomb_sec, aes(x=position, y=rate)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 19") +
  scale_x_continuous(breaks=c(109000000,109200000,109400000,109600000,109800000,110000000),
                     labels=c("109.0Mb", "109.2Mb", "109.4Mb", "109.6Mb", "109.8Mb", "110.0Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


#plot out association plots (2x), genes and recombinations rate
p1 <- (chr2_total/genes.plot/recomb_rate) + plot_layout(heights = c(3, 1, 1))

ggplot2::ggsave(
  "chr2_region_RA.jpg",
  width = 7,
  height = 5,
  dpi = 300
)



