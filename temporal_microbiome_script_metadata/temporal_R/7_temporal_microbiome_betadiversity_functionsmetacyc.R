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

##rarefied data- 2394389 (DRM5) ###################################################################################################################################################################################################################################################################################################
#Import data 
getwd() 

#import data
metacyc_table <- read.table("metacyc_rarefied_table_2394389.txt", sep = '\t', row.names = 1, header = T, strip.white = T) 
metacyc_table <- t(metacyc_table) # transpose to wide format
rowSums(metacyc_table) 

metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")

metacyc_table_metadata <- cbind(metacyc_table, metadata)
# Beta Diversity analysis##########################################################################################################################################################################################################################################################################################################################################################################################################################

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################

adonis.metacyc.bray <- adonis(metacyc_table_metadata[,1:411] ~ metacyc_table_metadata$day.time, method = "bray")
adonis.metacyc.bray$aov.tab

pairwiseAdonis::pairwise.adonis(metacyc_table_metadata[,1:411], metacyc_table_metadata$day.time, 
                                p.adjust.m ='fdr', sim.method = 'bray') 

#2. PERMDISP######################################################################################################################################################################################################################################################################################
rarefied_bray_dist <- vegdist(metacyc_table_metadata[,1:411], method = "bray") 
bdisp <- betadisper(rarefied_bray_dist,metacyc_table_metadata$day.time, type=c("centroid"))
bdisp  
aov.bdisp <-anova(bdisp)
aov.bdisp      
permutest(bdisp) 

permutest(bdisp, permutations = 999, pairwise = TRUE) 

labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp, cex=1, pch=15:17,
     main="Bray-Curtis Dissimilarity index", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)  


# PCoA plotting #########################################################################################################################################################################################################################################################################################
# calculate principal coordinates analysis (Bray-Curtis)
bray.pcoa.eigen <- cmdscale(rarefied_bray_dist, k = 2, eig = T)

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
bray.pcoa.eigen.plotting <- as.data.frame(bray.pcoa.eigen$points)
colnames(bray.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
bray.pcoa.eigen.plotting$sample <- rownames(bray.pcoa.eigen.plotting)

bray.pcoa.eigen.plotting <- cbind(bray.pcoa.eigen.plotting, metadata)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
(bray.pcoa.eigen$eig[1]/(sum(bray.pcoa.eigen$eig)))*100 #69.53%

(bray.pcoa.eigen$eig[2]/(sum(bray.pcoa.eigen$eig)))*100 #7.84%

# create a PCoA plot
pcoa.meio.bray.plot <- ggplot(bray.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = day.time)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (69.53%)") +
  ylab("PCoA 2 (7.84%)") +
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
pcoa.meio.bray.plot + annotate("text", x = -0.2, y = 0.08, label = "PERMANOVA:\np-value > 0.05\n PERMDISP:\np-value < 0.05")

#ggsave(filename = "metacyc_rarefied_bray_pcoa_plot", width = 20, height = 15, units = "cm", device = "pdf", dpi = "print")

