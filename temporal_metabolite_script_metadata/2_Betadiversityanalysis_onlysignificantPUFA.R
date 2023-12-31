getwd()

############Download gg3D###########################################################################################################################################################################################################################################################################################################
devtools::install_github("AckerDWM/gg3D")

#######################################################################################################################################################################################################################################################################################################################################
## Load packages
library(vegan)
library(cluster)
library(pairwiseAdonis)
library(dplyr)
library(ape)
library(ggplot2)
library(gg3D) #gg3D is a package created to extend ggplot2 to produce 3D plots

#Import data 
metabolites_all <- read.delim("metabolite_products_rawdata.txt", sep = '\t',row.names = 1 , header = T, strip.white = T) #this version follows the data structure in the file: "asv_table_wide_rarefied5612"
metabolites_all <- metabolites_all[-c(8,18),] #remove the BRM3 and DRM3 samples

#Use the significant 4 from the ANOVA test
metabolites_all_significantfour <- metabolites_all[, c("X5.F2t.IsoP", "X10.F4.NeuroP", "X13.F4.NeuroP", "X16.F1.PhytoP")]

metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")
metadata <- metadata[-c(8,18),] #remove BRM3 and DRM3 samples

##########################################################################################################################################################################################################################################################################################################################################################################################################################
## Statistical analysis on dist matrix ----
## Combine the data first ###
metabolites_four_metadata <- cbind(metabolites_all_significantfour, metadata)
metabolites_four_metadata$date <- factor(metabolites_four_metadata$date)
metabolites_four_metadata$day.time <- factor(metabolites_four_metadata$day.time )

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! it could be not. 2. The result always changes because of the random permutation and so the significance return will always be different. 
adonis.nonenzymaticmetabolites <- adonis(metabolites_four_metadata[,1:4] ~ metabolites_four_metadata$day.time, method = "euclidean")
adonis.nonenzymaticmetabolites$aov.tab

#                                                       Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#metabolites_all.no7HDHA.nonenzymatic_metadata$day.time  3   9226656 3075552  4.0186 0.46269  0.013 *
#Residuals                                              14  10714728  765338         0.53731         
#Total                                                  17  19941384                 1.00000       

#Post-hoc test: I plan to use fdr as correction method (used in QIIME2 as well)
pairwise.adonis(metabolites_four_metadata[,1:4], metabolites_four_metadata$day.time, 
                                p.adjust.m ='fdr', sim.method = 'euclidean') 

#2. PERMDISP######################################################################################################################################################################################################################################################################################
metabolites_four_euclidean_dist <- vegdist(metabolites_four_metadata[,1:4], method = "euclidean") 

bdisp <- betadisper(metabolites_four_euclidean_dist, metabolites_four_metadata$day.time, type=c("centroid"))
bdisp   
aov.bdisp <-anova(bdisp)
aov.bdisp       #aov and permutest (next code) produced the same result - with minor difference in the p-value -see Pat Schloss video on the interpretation CC208
permutest(bdisp, pairwise = TRUE) #Significant : How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#Response: Distances
#              Df   Sum Sq Mean Sq      F N.Perm Pr(>F)   
#Groups          3 1575609  525203 3.5135    999  0.005 **
#Residuals      14 2092729  149481  

#Simple overview with  PCoA plot
labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp, cex=1, pch=15:17,
     main="Euclidean Dissimilarity index", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)  #Even the confidence here chosen is 0.68, so I guess my code above for nmds uses 0.7 is alright because we mainly only want to visualised. -not with any statistics logic included for the ellipse


#PCoA plotting #########################################################################################################################################################################################################################################################################################
#Following the code from here: https://www.rpubs.com/roalle/mres_2019 

# calculate principal coordinates analysis (Euclidean)
euclidean.pcoa.eigen <- cmdscale(metabolites_four_euclidean_dist, k = 2, eig = T)

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
euclidean.pcoa.eigen.plotting <- as.data.frame(euclidean.pcoa.eigen$points)
colnames(euclidean.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
euclidean.pcoa.eigen.plotting$sample <- rownames(euclidean.pcoa.eigen.plotting)

euclidean.pcoa.eigen.plotting <- cbind(euclidean.pcoa.eigen.plotting, metadata)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
(euclidean.pcoa.eigen$eig[1]/(sum(euclidean.pcoa.eigen$eig)))*100 #91.83845%

(euclidean.pcoa.eigen$eig[2]/(sum(euclidean.pcoa.eigen$eig)))*100 #4.355621%

# create a PCoA plot - result same as QIIME2!
pcoa.meio.bray.plot <- ggplot(euclidean.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = day.time)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (91.83%)") +
  ylab("PCoA 2 (4.36%)") +
  labs(color = "Timepoints") +
  scale_color_hue(labels = c("T_0915", "T_1230", "T_1612", "T_1740")) +
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), #changing the title of legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) 
        
pcoa.meio.bray.plot + annotate("text", x = -2000, y = -800, label = "PERMANOVA: p-value < 0.05\n PERMDISP: p-value < 0.05")
#ggsave(filename = "nonenzymaticmetabolites_pcoa_foursignificantmetabolites", width = 20, height = 15, units = "cm", device = "pdf", dpi = "print")
