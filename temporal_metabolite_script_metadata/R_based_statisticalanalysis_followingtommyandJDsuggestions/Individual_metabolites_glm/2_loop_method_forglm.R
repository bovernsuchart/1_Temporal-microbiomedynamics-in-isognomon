# glm for those test that violated either one of the assumptions - with loop method

library(ggplot2)
library(car)
library(dplyr)
library(DHARMa) #same as used by JD
library(emmeans)
library(report)
library(multcomp)

#Import data metabolites data###################################################################################################################################
metadata <- read.delim("metabolite_products_rawdata_forindividualanalsysis_inR.txt", header = T, row.names = 1,comment.char="#")
metadata$day.time <- factor(metadata$day.time)
colnames(metadata) #do this to get an accurate name of the metabolites in the table

#Select the groups that needs to do glm
colnames(metadata)
glm_metabolite <- c("day.time", "X17.Dihomo.IsoP", "X13.F4.NeuroP", "X9.F1.PhytoP", "X16.B1.PhytoP", "X15.F3t.IsoP")

glm_metabolite_metadata <- metadata[colnames(metadata) %in% glm_metabolite]

## Individual metabolites analysis with glm using gamma(link = "inverse') #############################################################################################################################
#loop method
glm_model <- list()

for ( i in 2:6) { 
  glm_model[[i]] <- glm(glm_metabolite_metadata [[i]] ~ day.time, family = Gamma(link = "inverse"), glm_metabolite_metadata )
    #glm_modelsimulation[[i]] <- simulateResiduals(glm_model[[i]], n=250, plot=T) #for residual diagnostics
  }

#residuals diagnostics - very slow - better to do one by one for this
#glm_modelsimulation <- list()
 
#for ( x in 2:6) { 
#glm_modelsimulation[[x]] <- simulateResiduals(glm_model[[x]], n=250, plot=T) #for residual diagnostics
#}


#Data presentation#####
#Adrian method - using Anova() - which referenced from here: Fox & Weisberg (2019) An R companion to applied regression. See the paper by Eckert et al. (2021) and Costantini et al. (2021) - "Wild gut microbiomes reveal individuals, species, and location as drivers of variation in two critically endangered Hawaiian honeycreepers"
#?Anova()
Anova_model <- list()
Anova_model_pvalue <- vector() #to extract the pvalue
summary_glmmodel <- list()
report_glmmodel <- list()

for (j in 2:6){ 
  Anova_model[[j]] <- Anova(glm_model[[j]], test.statistic = "LR") #preferred
  Anova_model_pvalue[[j]] <- Anova_model[[j]]$`Pr(>Chisq)`
  summary_glmmodel[[j]] <- summary(glm_model[[j]]) #JD method
  report_glmmodel[[j]] <- report(glm_model[[j]]) #JD method
}


#Extract the pvalue model
LRChiseq <- vector()
Df <- vector()
Chisq_pvalue <- vector()

for (y  in 2: 6) {
   LRChiseq[[y]] <-  Anova_model[[y]]$`LR Chisq`
   Df[[y]] <-  Anova_model[[y]]$Df
   Chisq_pvalue[[y]] <- Anova_model[[y]]$`Pr(>Chisq)`
}

#compile all the data together with metabolites data#################################################3
glm_metabolite <- c("day.time", "X17.Dihomo.IsoP", "X13.F4.NeuroP", "X9.F1.PhytoP", "X16.B1.PhytoP", "X15.F3t.IsoP")

glm_result <- data.frame(cbind(glm_metabolite, LRChiseq, Df, Chisq_pvalue))



#Pairwise comparison## _ Adrian method better! Which is shown in some papers! Agostini et al. (2021) & Wood et al. (2022)################################
emmeans <- list()
glht <- list()

anova_deviance <- list()
for (k in 2:6){
  #emmeans result will be the same as glht
  emmeans[[k]] <- emmeans(glm_model[[k]], pairwise ~ day.time  ) # Pairwise comparison - Adrian added - better to use this t-ratio since my data set is small
  glht[[k]] <- summary(glht(glm_model[[k]], mcp(day.time="Tukey"))) #Follow Eckert et al. (2021)
  anova_deviance[[k]] <- anova(glm_model[[k]], test = "Chisq")  # (Extra) JD anova method - uses analysis of deviance table
  }


