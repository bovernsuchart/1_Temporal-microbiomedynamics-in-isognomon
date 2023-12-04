# glm for those test that violated either one of the assumptions

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

## Individual metabolites analysis with glm using gamma(link = "inverse') #############################################################################################################################

#4. X17.Dihomo.IsoP ################################################################################################################################ X17.Dihomo.IsoP
glm4 <- glm(X17.Dihomo.IsoP ~ day.time, family = Gamma(link="inverse"), metadata)
summary(glm4)
report(glm4)

# Residual diagnostic
simulateResiduals(glm4, n=250, plot=T)

#Adrian method - using Anova() - which referenced from here: Fox & Weisberg (2019) An R companion to applied regression. See the paper by Eckert et al. (2021) and Costantini et al. (2021) - "Wild gut microbiomes reveal individuals, species, and location as drivers of variation in two critically endangered Hawaiian honeycreepers"
#?Anova()
Anova(glm4, test.statistic = "LR")

# Pairwise comparison - Adrian added - better to use this t-ratio since my data set is small
emmeans(glm4, pairwise ~ day.time) #The result is the same as glht

#Follow Eckert et al. (2021)
summary(glht(glm4, mcp(day.time="Tukey"))) #The result is the same as emmeans

# (Extra) JD anova method - uses analysis of deviance table
anova(glm4, test = "Chisq") 

#7. X13.F4.NeuroP ###############################################################################################################################
glm7 <- glm(X13.F4.NeuroP ~ day.time, family = Gamma(link="inverse"), metadata)
summary(glm7)
report(glm7)

# Residual diagnostic
simulateResiduals(glm7, n=250, plot=T)

#Adrian method - using Anova() - which referenced from here: Fox & Weisberg (2019) An R companion to applied regression. See the paper by Eckert et al. (2021) and Costantini et al. (2021) - "Wild gut microbiomes reveal individuals, species, and location as drivers of variation in two critically endangered Hawaiian honeycreepers"
#?Anova()
Anova(glm7, test.statistic = "LR")

# Pairwise comparison - Adrian added - better to use this t-ratio since my data set is small
emmeans(glm7, pairwise ~ day.time) #The result is the same as glht

#Follow Eckert et al. (2021)
summary(glht(glm7, mcp(day.time="Tukey"))) #The result is the same as emmeans

# (Extra) JD anova method - uses analysis of deviance table
anova(glm7, test = "Chisq") 

#10. X9.F1.PhytoP ###############################################################################################################################
glm10 <- glm(X9.F1.PhytoP ~ day.time, family = Gamma(link="inverse"), metadata)
summary(glm10)
report(glm10)

# Residual diagnostic
simulateResiduals(glm10, n=250, plot=T)

#Adrian method - using Anova() - which referenced from here: Fox & Weisberg (2019) An R companion to applied regression. See the paper by Eckert et al. (2021) and Costantini et al. (2021) - "Wild gut microbiomes reveal individuals, species, and location as drivers of variation in two critically endangered Hawaiian honeycreepers"
#?Anova()
Anova(glm10, test.statistic = "LR")

# Pairwise comparison - Adrian added - better to use this t-ratio since my data set is small
emmeans(glm10, pairwise ~ day.time) #The result is the same as glht

#Follow Eckert et al. (2021)
summary(glht(glm10, mcp(day.time="Tukey"))) #The result is the same as emmeans

# (Extra) JD anova method - uses analysis of deviance table
anova(glm10, test = "Chisq") 

#12. X16.B1.PhytoP ###############################################################################################################################
glm12 <- glm(X16.B1.PhytoP ~ day.time, family = Gamma(link="inverse"), metadata)
summary(glm12)
report(glm12)

# Residual diagnostic
simulateResiduals(glm12, n=250, plot=T)

#Adrian method - using Anova() - which referenced from here: Fox & Weisberg (2019) An R companion to applied regression. See the paper by Eckert et al. (2021) and Costantini et al. (2021) - "Wild gut microbiomes reveal individuals, species, and location as drivers of variation in two critically endangered Hawaiian honeycreepers"
#?Anova()
Anova(glm12, test.statistic = "LR")

# Pairwise comparison - Adrian added - better to use this t-ratio since my data set is small
emmeans(glm12, pairwise ~ day.time) #The result is the same as glht

#Follow Eckert et al. (2021)
summary(glht(glm12, mcp(day.time="Tukey"))) #The result is the same as emmeans

# (Extra) JD anova method - uses analysis of deviance table
anova(glm12, test = "Chisq") 


#14. X15.F3t.IsoP ###############################################################################################################################
glm14 <- glm(X15.F3t.IsoP ~ day.time, family = Gamma(link = "inverse"), metadata)
summary(glm14)
report(glm14)

# Residual diagnostic
simulateResiduals(glm14, n=250, plot=T)

#Adrian method - using Anova() - which referenced from here: Fox & Weisberg (2019) An R companion to applied regression. See the paper by Eckert et al. (2021) and Costantini et al. (2021) - "Wild gut microbiomes reveal individuals, species, and location as drivers of variation in two critically endangered Hawaiian honeycreepers"
#?Anova()
Anova(glm14, test.statistic = "LR")

# Pairwise comparison - Adrian added - better to use this t-ratio since my data set is small
emmeans(glm14, pairwise ~ day.time) #The result is the same as glht

#Follow Eckert et al. (2021)
summary(glht(glm14, mcp(day.time="Tukey"))) #The result is the same as emmeans

# (Extra) JD anova method - uses analysis of deviance table
anova(glm14, test = "Chisq") 

