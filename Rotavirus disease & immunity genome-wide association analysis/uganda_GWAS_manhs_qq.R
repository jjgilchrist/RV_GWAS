rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)

#read in summary statistics: GWAS of anti-RV IgA levels in 2 and 3 year olds in Uganda
#summary statistics are randomly thinned (~1 in 10) with the exception of chromosomes and 19
#full summary statistics will be deposited with the GWAS Catalog (https://www.ebi.ac.uk/gwas/)

disc.add <- read.table("uganda_GWAS_summ_stats.thinned.txt.gz", header = T)

disc.add$chr <- factor(disc.add$chr, levels = c(1:22))
levels(disc.add$chr) <- c(1:22)

#need to highlight secretor and chr2 locus
 
disc.add$label <- NA
disc.add$label[which(disc.add$rsid=="19:49206674:G:A")] <- "rs601338"
disc.add$label[which(disc.add$rsid=="2:109558305:C:T")] <- "rs59241810"

disc.add$anno <- NA
disc.add$anno[which(disc.add$rsid=="19:49206674:G:A")] <- 1
disc.add$anno[which(disc.add$rsid=="2:109558305:C:T")] <- 1


#cumulative base pair position

disc.add2 <- disc.add %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(disc.add, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = disc.add2 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#set colour palette
cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")

#plot manhattan plot
uganda.manh <- ggplot(disc.add2, aes(x=BPcum, y=-log10(p))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  
  
  # Add test using ggrepel to avoid overlapping
  geom_text_repel( data=subset(disc.add2, anno==1), aes(label=label), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 2) +
  
  #gwas sig line
  geom_segment(aes(x = 10177, y = 7.30103, xend = 2879943885, yend = 7.30103), data = disc.add2, linetype="dashed", color = "red") +


  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","2","4","6","8","10","12","14","16","18"), breaks = c(0,2,4,6,8,10,12,14,16,18), expand = c(0, 0), limits= c(0,18)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20)
  )

ggsave(
  "uganda.manh.jpg",
  width = 14,
  height = 3.5,
  dpi = 300
)

#compile observed and expected (null sumamry statistics)

df <- data.frame(observed = -log10(sort(na.omit(disc.add$p))),
    expected = -log10(ppoints(length(na.omit(disc.add$p)))))
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
#calculate genomic inflation/lambda
  chisq <- qchisq(1-disc.add$p,1)
  round(median(chisq)/qchisq(0.5,1),2)

  lambda <- paste("lambda == ", round(median(chisq)/qchisq(0.5,1),2))

#plot QQ plot
rv.qq <- ggplot(df) +
  geom_abline(intercept = 0, slope = 1, size = 1, colour = cols3[6]) +
  geom_point(aes(expected, observed), size = 3, colour = cols3[2]) +
  xlab(log10Pe) +
  ylab(log10Po) +
  annotate("text", x = 1, y = 16, label = lambda, parse = TRUE, size = 10) +
  theme_bw() + ylim(NA, 18) + 
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20), plot.title=element_text(size=30, face = "bold", hjust = 0.5))

ggsave(
  "uganda.qq.jpg",
  width = 7,
  height = 7,
  dpi = 300
)

