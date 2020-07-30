library(ggplot2)
library("dplyr")
library("tidyr")
library(patchwork)
library("igraph")
library(calibrate)
library("umap")

library(viridis)

#########################################################################
### REACTOME & viridis color

setwd("c://Users/egamazon/Dropbox/Papers/GTEx_v8_Network_Cambridge/Supplementary/")
a = read.table('data/Hippocampus_community.txt', header=T, sep="\t")

ggplot(data.frame(a), aes(x=reorder(gsub("Homo.*","",Term), -P.value), y=-log10(P.value)))  +
theme(panel.background = element_rect(fill="grey", color="black"), panel.grid.major = element_line(color="grey90"), panel.grid.minor = element_line(color="grey90")) +
theme(axis.text=element_text(size=22), axis.title=element_text(size=30), text=element_text(size=22))   + xlab(expression("Reactome")) +
  geom_bar(stat="identity", fill="darkgreen") + ylab(expression("-log(P-value)")) + 
  geom_col(aes(fill = Adjusted_P.value)) + 
	scale_fill_viridis() +
coord_flip() + geom_hline(yintercept=-log10(0.05), color="red") +guides(fill = guide_legend(reverse = FALSE)) + labs(fill="Adjusted P-value") 

#########################################################################
### Correlation plot
b = read.table('data/15_member_hippocampal_community.v1.txt', header=T)
rownames(b) = b[,1]

b[is.na(b)] <- 0

c = ( b[,-1] + t(b[,-1]) ) 
c[c >1.98] = Inf

c$y = rownames(c)
c <- c %>% gather(x, value, colnames(c)[1:15])

# centered vjust = 0.5
ggplot(c, aes(x, y, fill=value)) + theme(axis.text=element_text(size=20), axis.title=element_text(size=30), text=element_text(size=22))  + 
   xlab(expression(bold(""))) + ylab(expression(bold(""))) +
  geom_tile() + # scale_fill_gradient(low = "white", high = "orange") +  
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + # labs(tag = "B")
  scale_fill_viridis() 


#########################################################################
### Network

c = ( b[,-1] + t(b[,-1]) ) 
c[c >1.98] = 1

mat = data.matrix(c) 

mat[mat < 0.8] <- 0

network = graph_from_adjacency_matrix( mat, mode="undirected", diag=F, weighted=T)

plot(network, vertex.size=30,  vertex.label.font=1,   vertex.color="yellow" ,layout=layout_with_kk,
 vertex.label.color = "darkgreen", vertex.label.cex=1.3, vertex.label.pos=5, vertex.shape="circle", edge.width=6)


#########################################################################
### Stats on communities

a = read.table('data/CommunitySize.txt', header=T, sep="\t")

wilcox.test(a[c(7:19),]$std, a[-c(7:19),]$std)

par(mar=c(5,6,4,1)+.1)

par(mfrow=c(1,2))
hist(a$n_comm, main="", cex=2.0, cex.lab=2.0, cex.axis=2.0, xlab="Number of Communities", col="darkgreen", xlim=c(0, 260), breaks=20)
#plot(a$mean, a$std, main="", cex=2.0, cex.lab=2.0, cex.axis=2.0, xlab="Mean", ylab="Standard Deviation", col="red")
plot(a$n_comm, a$mean, main="", cex=2.0, cex.lab=2.0, cex.axis=2.0, xlab="Number of Communities", ylab="Mean Community Size", col="red", ylim=c(0, 9), xlim=c(25,300))
textxy(a$n_comm, a$mean, labs=ifelse((a$n_comm > 150) | (a$mean > 7.25), as.character(a$Tissue), ""), col="blue", cex=1.2, pos=3)  

###########################################################################
### Summary data on communities

d = read.table('data/CommunitySize.txt', header=T, sep="\t")

p1 = ggplot(data.frame(d), aes(n_comm)) + 
    geom_histogram(data = data.frame(d),  binwidth=10, fill="darkgreen", col="black") + 
    xlab(expression("Number of Communities")) + ylab(expression("Frequency")) +
    theme(plot.title = element_text(hjust = 0.5),  axis.title=element_text(size=30), text=element_text(size=30), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) + labs(tag = expression("a")) + theme( axis.line = element_line(colour = "grey90", size = 2, linetype = "solid")) +
	xlim(c(25, 300)) 

p2 = ggplot(data.frame(d), aes(x=n_comm, y=mean)) + geom_point(size=6, col="darkgreen") + geom_text(aes(label=ifelse((n_comm > 150) , as.character(Tissue), "")),hjust=1, vjust=1, size=8) + xlab(expression("Number of Communities")) + ylab(expression("Mean Size")) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.title=element_text(size=30), text=element_text(size=30), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) +
 xlim(c(25, 300)) + ylim(c(0, 8)) + labs(tag=expression("b")) + theme( axis.line = element_line(colour = "grey90", size = 2, linetype = "solid")) 

# p1+p2

p3 = ggplot(d, aes(x = Tissue, y = 18365 - ngenes_no_comm)) + geom_bar(stat="identity", aes(fill=  18365 - ngenes_no_comm)) + scale_fill_viridis() + theme_minimal() +  labs(tag=expression("c")) +
    theme(panel.border = element_blank(), legend.position="none", plot.title = element_text(hjust = 0.5), axis.title=element_text(size=30), text=element_text(size=26), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme( axis.line = element_line(colour = "grey90", size = 2, linetype = "solid")) +
  xlab(expression("Tissue")) + ylab(expression("Frequency")) +  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 


(p1 + p2) / p3

###########################################################################
####  Variational Autoencoder results

clin.tcga = read.table('data/clinical_data.TCGA.tsv', header=T, fill=T, sep="\t")
enc  = read.table('data/encoded_rnaseq_onehidden_warmup_batchnorm.tsv', fill=T, sep="\t")
comb=merge(clin.tcga, enc, by.x=1, by.y=2)

####  Metastasis analysis

comb1 = subset(comb, comb[,9] %in% c("Metastatic", "Additional Metastatic", "Primary Tumor", "Solid Tissue Normal"))

vaeplot<-ggplot(data=comb1, aes(x=sample_type, y=V6)) + geom_boxplot(aes(fill="darkgreen")) + geom_point() + theme_minimal() +  labs(tag=expression(bold("C"))) +
    theme(panel.border = element_blank(), legend.position="none", plot.title = element_text(hjust = 0.5), axis.title=element_text(size=22), text=element_text(size=18), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme( axis.line = element_line(colour = "grey90", size = 2, linetype = "solid")) +
 theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab("Sample Type") + ylab("Latent Feature 4 Activation") + 
	scale_fill_brewer(palette="Dark2") # + geom_jitter(shape=16, position=position_jitter(0.2))

vaeplot 

### Compare
wilcox.test(comb1[comb1$sample_type=="Primary Tumor",]$V6, comb1[comb1$sample_type=="Metastatic",]$V6)

# data:  comb1[comb1$sample_type == "Primary Tumor", ]$V6 and comb1[comb1$sample_type == "Metastatic", ]$V6
# W = 632560, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(comb1[comb1$sample_type=="Solid Tissue Normal",]$V6, comb1[comb1$sample_type=="Metastatic",]$V6)

# data:  comb1[comb1$sample_type == "Solid Tissue Normal", ]$V6 and comb1[comb1$sample_type == "Metastatic", ]$V6
# W = 64866, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

###########################################################################
### Stemness

stem = read.table('data/VAE_with_StemnessScore.v1.txt', header=T, sep="\t")
 
stemness =  ggplot(stem, aes(x=log2(X4), y=log2(DNAss) )) + geom_point(size=6, col="yellow")  + xlab(expression("Latent Feature 4 Activation")) + ylab(expression("Stemness Index")) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.title=element_text(size=22), text=element_text(size=18), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) +
  labs(tag=expression(bold("D"))) + theme( axis.line = element_line(colour = "grey90", size = 2, linetype = "solid")) +
#  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, col="darkgreen")
