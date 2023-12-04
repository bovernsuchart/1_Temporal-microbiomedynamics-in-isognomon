getwd()

############Download gg3D#############################################################################################################################################################################################################################################################################################################
devtools::install_github("AckerDWM/gg3D")

######################################################################################################################################################################################################################################################################################################################################
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

metadata <- read.delim("metadata_daynight.txt", header = T, row.names = 1,comment.char="#")
metadata <- metadata[-c(8,18),] #remove BRM3 and DRM3 samples

##########################################################################################################################################################################################################################################################################################################################################################################################################################
## Statistical analysis on dist matrix ----
## Combine the data first ###
metabolites_all_metadata <- cbind(metabolites_all, metadata)

metabolites_all_metadata$date <- factor(metabolites_all_metadata$date)
metabolites_all_metadata$day.time <- factor(metabolites_all_metadata$day.time )

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! it could be not. 2. The result always changes because of the random permutation and so the significance return will always be different. 
adonis.nonenzymaticmetabolites <- adonis(metabolites_all_metadata[,1:15] ~ metabolites_all_metadata$day.time, method = "euclidean")
adonis.nonenzymaticmetabolites$aov.tab

#                                                       Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#metabolites_all.no7HDHA.nonenzymatic_metadata$day.time  3  42242237 14080746  2.1423 0.31463  0.086 .
#Residuals                                              14  92019170  6572798         0.68537         
#Total                                                  17 134261407                  1.00000    

#Post-hoc test: I plan to use fdr as correction method (used in QIIME2 as well)
pairwise.adonis(metabolites_all_metadata[,1:15], metabolites_all_metadata$day.time, 
                                p.adjust.m ='fdr', sim.method = 'euclidean') 

#2. PERMDISP######################################################################################################################################################################################################################################################################################
metabolites_all_metadata_euclidean_dist <- vegdist(metabolites_all_metadata[,1:15], method = "euclidean") 

bdisp <- betadisper(metabolites_all_metadata_euclidean_dist, metabolites_all_metadata$day.time, type=c("centroid"))
bdisp   
aov.bdisp <-anova(bdisp)
aov.bdisp       #aov and permutest (next code) produced the same result - with minor difference in the p-value -see Pat Schloss video on the interpretation CC208
permutest(bdisp, pairwise = TRUE) #Significant : How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#Response: Distances
#              Df   Sum Sq Mean Sq      F N.Perm Pr(>F)   
#Groups        3 14866469 4955490 6.8803    999  0.004 **
#Residuals    14 10083393  720242   

#Simple overview with  PCoA plot
labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp, cex=1, pch=15:17,
     main="Euclidean Dissimilarity index", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)  


#PCoA plotting #########################################################################################################################################################################################################################################################################################
#Following the code from here: https://www.rpubs.com/roalle/mres_2019 

# calculate principal coordinates analysis (Euclidean)
euclidean.pcoa.eigen <- cmdscale(metabolites_all_metadata_euclidean_dist, k = 2, eig = T)

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
euclidean.pcoa.eigen.plotting <- as.data.frame(euclidean.pcoa.eigen$points)
colnames(euclidean.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
euclidean.pcoa.eigen.plotting$sample <- rownames(euclidean.pcoa.eigen.plotting)

euclidean.pcoa.eigen.plotting <- cbind(euclidean.pcoa.eigen.plotting, metadata)

#Calculate the proportion of variance in the data which is explained by the first two PCoA axes
(euclidean.pcoa.eigen$eig[1]/(sum(euclidean.pcoa.eigen$eig)))*100 #78.79866%

(euclidean.pcoa.eigen$eig[2]/(sum(euclidean.pcoa.eigen$eig)))*100 #14.75%

# create a PCoA plot - result same as in QIIME2!
pcoa.meio.bray.plot <- ggplot(euclidean.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = day.time)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (78.80%)") +
  ylab("PCoA 2 (14.75%)") +
  labs(color = "Timepoints") +
  scale_color_hue(labels = c("T_0915", "T_1230", 'T_1612', "T_1740")) +
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), #changing the title of legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) 
        
pcoa.meio.bray.plot + annotate("text", x = 9750, y = 3750, label = "PERMANOVA: p-value > 0.05\n PERMDISP: p-value < 0.05")
#ggsave(filename = "nonenzymaticmetabolites_pcoa_allmetabolites", width = 20, height = 15, units = "cm", device = "pdf", dpi = "print")

