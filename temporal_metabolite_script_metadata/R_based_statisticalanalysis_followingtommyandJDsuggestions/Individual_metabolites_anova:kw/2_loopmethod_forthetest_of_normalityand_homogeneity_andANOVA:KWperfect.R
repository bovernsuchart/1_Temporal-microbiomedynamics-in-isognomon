#Individual metabolites analysis in R (using loop function) - 13/02/2023

## Load packages
library(ggplot2)
library(FSA)
library(phyloseq)
library(forcats)
library(reshape2)
library(dplyr) 
library(carData)
library(car)
library(multcomp)
library(permute)
library(lattice)
library(vegan)
library(performance)
library(cowplot)
library(vegan)

#Import data metabolites data###################################################################################################################################
metadata <- read.delim("metabolite_products_rawdata_forindividualanalsysis_inR.txt", header = T, row.names = 1,comment.char="#")
metadata$day.time <- factor(metadata$day.time)
colnames(metadata) #do this to get an accurate name of the metabolites in the table

## Individual metabolites analysis #############################################################################################################################

#Linear model, Normality and Homogeneity of Variance
model <- list() #need to specify the data type for the object first
normality <- list() 
homo.variance <- list()

for (i in 2:16) {
    model[[i]] <- lm(metadata[[i]] ~ day.time, metadata)
    normality[[i]] <- shapiro.test(residuals(object = model[[i]]))
    homo.variance[[i]] <- leveneTest(model[[i]], center = mean)
}

#Extract normality's and homogeneity's p-value
normality_pvalue <- vector()
homo.variance_pvalue <- vector()
for (y in 2:16){
  normality_pvalue[[y]] <- normality[[y]]$p.value
  homo.variance_pvalue[[y]] <- homo.variance[[y]][1,]$`Pr(>F)` #this is a data frame so the code is like this
  }

#lm model in ANOVA table 
anova_table <- list()
for (j in 2:16){
   anova_table[[j]] <- anova(lm(metadata[[j]] ~ day.time, metadata))
  }

#Extract model p-value
model_pvalue <- vector()
for (x in 2:16){
  model_pvalue[[x]] <- anova_table[[x]][1,]$`Pr(>F)`  #or you can use this code here: anova_table[[x]][1,5]
}

#KW test & extract the p-value
KW <- list()
KW_pvalue <- vector()
for (k in 2:16){
  KW[[k]] <- kruskal.test(metadata[[k]] ~ day.time, metadata)
  KW_pvalue[[k]] <- KW[[k]]$p.value
}


#remove all the NA values and create new object
model_pvalue #with NA value
model_pvalue2 <- round(model_pvalue[!is.na(model_pvalue)], digits = 4) #remove the NA value 
normality_pvalue2 <- round(normality_pvalue[!is.na(normality_pvalue)], digits = 4)
homo.variance_pvalue2 <- round(homo.variance_pvalue[!is.na(homo.variance_pvalue)], digits = 4)
KW_pvalue2 <- round(KW_pvalue[!is.na(KW_pvalue)], digits = 4) #remove the NA value 


#compile all the data together with metabolites data
metabolites <- colnames(metadata[2:16])
metabolites
model_information <- data.frame(cbind(metabolites, model_pvalue2, normality_pvalue2, homo.variance_pvalue2, KW_pvalue2))

