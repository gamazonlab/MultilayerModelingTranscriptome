library(ggplot2)
library(dplyr)
library(patchwork)

gwas.dat = read.table('data/CRP_st_Whole_Blood.txt', header=T)
gwas.dat = subset(gwas.dat, !is.na(gwas.dat$CHR))

gwas.dat = gwas.dat[order(gwas.dat$CHR),]

cgenes = read.table('data/community_genes.txt', header=F)

gwas.dat = subset(gwas.dat, gwas.dat$gene_name %in% cgenes[,1])


nCHR <- length(unique(gwas.dat$CHR))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(gwas.dat$CHR)){
  nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$POS)
  gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"POS"] + s
  s <- s + nbp[i]
}

#
axis.set <- gwas.dat %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(gwas.dat$pvalue)))) + 2 


sig = 0.05/ dim(gwas.dat)[1]

# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

manhattan <- ggplot(gwas.dat, aes(x = BPcum, y = -log10(pvalue), 
                                 color = as.factor(CHR), size = -log10(pvalue))) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
#  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_color_manual(values = rep(c("red", "blue"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = expression(paste("-log"[10], plain(P)))) + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 18, vjust = 0.5),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 22)
  ) +  labs(tag = expression(bold("B"))) + geom_text(aes(label=ifelse(pvalue < 1e-10, as.character(gene_name),"")), size=5, check_overlap=T, nudge_x = 1.25, nudge_y = -2.25)

 manhattan
