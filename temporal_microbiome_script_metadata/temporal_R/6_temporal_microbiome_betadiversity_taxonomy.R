## Load packages

library(permute)
library(dplyr)
library(lattice)
library(vegan)
library(cluster)
library(pairwiseAdonis)
library(ape)
library(ggplot2)
library(gg3D) 

##rarefied data- 21282###################################################################################################################################################################################################################################################################################################
#Import data 

getwd() 
asv_table_rarefied21282 <- read.delim("rarefied_feature_table-21282.txt", sep = '\t', header = T, row.names = 1, strip.white = T)
asv_table_wide_rarefied21282 <- t(asv_table_rarefied21282) # transpose to wide format

metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")

##########################################################################################################################################################################################################################################################################################################################################################################################################################
## Beta diversity analysis ----
#PERMANOVA 

#combine both the datasets
asv_rarefied21282_metadata <- cbind(asv_table_wide_rarefied21282, metadata)

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
adonis.microbiome <- adonis(asv_rarefied21282_metadata[,1:4586] ~ asv_rarefied21282_metadata$day.time, method = "bray")
adonis.microbiome$aov.tab

#Post-hoc test: 
pairwiseAdonis::pairwise.adonis(asv_rarefied21282_metadata[,1:4586], asv_rarefied21282_metadata$day.time, 
                                p.adjust.m ='fdr', sim.method = 'bray') 

#2. PERMDISP######################################################################################################################################################################################################################################################################################
rarefied21282_bray_dist <- vegdist(asv_rarefied21282_metadata[,1:4586], method = "bray") 
bdisp <- betadisper(rarefied21282_bray_dist, asv_rarefied21282_metadata$day.time, type=c("centroid"))
bdisp   
aov.bdisp <-anova(bdisp)
aov.bdisp       

permutest(bdisp, permutations = 999, pairwise = TRUE) #Significant : How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#Simple overview of samples clusters with ordination plot
labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp, cex=1, pch=15:17,
     main="Bray-Curtis Dissimilarity index", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)  #Even the confidence here chosen is 0.68, so I guess my code above for nmds uses 0.7 is alright because we mainly only want to visualised. -not with any statistics logic included for the ellipse


#PCoA plotting #########################################################################################################################################################################################################################################################################################

# calculate principal coordinates analysis (Bray-Curtis)
bray.pcoa.eigen <- cmdscale(rarefied21282_bray_dist, k = 2, eig = T)

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
bray.pcoa.eigen.plotting <- as.data.frame(bray.pcoa.eigen$points)
colnames(bray.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
bray.pcoa.eigen.plotting$sample <- rownames(bray.pcoa.eigen.plotting)

bray.pcoa.eigen.plotting <- cbind(bray.pcoa.eigen.plotting, metadata)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
(bray.pcoa.eigen$eig[1]/(sum(bray.pcoa.eigen$eig)))*100 #20.22653%

(bray.pcoa.eigen$eig[2]/(sum(bray.pcoa.eigen$eig)))*100 #15.93464%

# create a PCoA plot
pcoa.meio.bray.plot <- ggplot(bray.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = day.time)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (20.23%)") +
  ylab("PCoA 2 (15.93%)") +
  labs(color = "Timepoints") +
  scale_color_hue(labels = c("T_0915", "T_1230", 'T_1612', "T_1740")) +
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())
pcoa.meio.bray.plot + annotate("text", x = -0.6, y = 0.6, label = "PERMANOVA:\n p-value > 0.05\n PERMDISP:\n p-value < 0.05")

#ggsave(filename = "rarefied21282_bray_pcoa_plot-finalised", width = 20, height = 15, units = "cm", device = "pdf", dpi = "print")
