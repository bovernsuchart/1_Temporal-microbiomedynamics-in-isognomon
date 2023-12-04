#Individual metabolites analysis in R - 06/02/2023

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

#1. X5.F2t.IsoP ###############################################################################################################################
#Calculate statistics By day.time
model1 <- lm( X5.F2t.IsoP ~ day.time, metadata) #create a model
anova(model1) #Result: F-value(3,14): 5.0745 ; p-value: 0.01385 *

#Normality test -
par(mfrow=c(2,2))
plot(model1) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals1 <- residuals(object = model1)
shapiro.test(x = model.residuals1) #result: W = 0.92573, p-value = 0.1633

#Homogeneity of Variance - raw data
model.levene1 <- leveneTest(model1, center = mean)
model.levene1  #Result: p-value = 0.2725

#2. X15.F2t.IsoP ###############################################################################################################################
#Calculate statistics By day.time
model2 <- lm( X15.F2t.IsoP ~ day.time, metadata) #create a model
anova(model2) #Result: F-value(3,14): 1.715 ; p-value:  0.2096 *

#Normality test -
par(mfrow=c(2,2))
plot(model2) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals2 <- residuals(object = model2)
shapiro.test(x = model.residuals2) #result: W = 0.90324, p-value = 0.06546

#Homogeneity of Variance - raw data
model.levene2 <- leveneTest(model2, center = mean)
model.levene2  #Result: p-value = 0.05589

#kruskal-wallis test #In my thesis record, I used this test! I may need to change it to later
model.kruskal2 <-kruskal.test(X15.F2t.IsoP ~ day.time, data = metadata)
model.kruskal2  #Kruskal-Wallis chi-squared = 4.9298, df = 3, p-value = 0.177 
dunnTest(X15.F2t.IsoP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#3. X7.Dihomo.IsoP ###############################################################################################################################
#Calculate statistics By day.time
model3 <- lm( X7.Dihomo.IsoP ~ day.time, metadata) #create a model
anova(model3) #Result: F-value(3,14): 3.0994 ; p-value:  0.06108

#Normality test -
par(mfrow=c(2,2))
plot(model3) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals3 <- residuals(object = model3)
shapiro.test(x = model.residuals3) #result: W = 0.90324, p-value = 0.3632

#Homogeneity of Variance - raw data
model.levene3 <- leveneTest(model3, center = mean)
model.levene3  #Result: p-value =0.01234

#kruskal-wallis test #In my thesis record, I used this test! I may need to change it to later
model.kruskal3 <-kruskal.test(X7.Dihomo.IsoP  ~ day.time, data = metadata)
model.kruskal3  #Kruskal-Wallis chi-squared = 8.2719, df = 3, p-value = 0.04071
dunnTest(X7.Dihomo.IsoP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#4. X17.Dihomo.IsoP ###############################################################################################################################
#Calculate statistics By day.time
model4 <- lm(X17.Dihomo.IsoP ~ day.time, metadata) #create a model
anova(model4) #Result: F-value(3,14): 0.2122 ; p-value:  0.8863

#Normality test -
par(mfrow=c(2,2))
plot(model4) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals4 <- residuals(object = model4)
shapiro.test(x = model.residuals4) #result: W = 0.90324, p-value = 0.03467

#Homogeneity of Variance - raw data
model.levene4 <- leveneTest(model4)
model.levene4  #Result: p-value = 0.8185

#kruskal-wallis test #In my thesis record, I used ANOVA test! I may need to change it to later
model.kruskal4 <-kruskal.test(X17.Dihomo.IsoP ~ day.time, data = metadata)
model.kruskal4  #Kruskal-Wallis chi-squared = 1.4561, df = 3, p-value = 0.6924
dunnTest(X17.Dihomo.IsoP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#5. X4.F4.NeuroP ###############################################################################################################################
#Calculate statistics By day.time
model5 <- lm( X4.F4.NeuroP ~ day.time, metadata) #create a model
anova(model5) #Result: F-value(3,14): 3.1303 ; p-value:  0.05953

#Normality test -
par(mfrow=c(2,2))
plot(model5) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals5 <- residuals(object = model5)
shapiro.test(x = model.residuals5) #result: W = 0.91931, p-value = 0.1257

#Homogeneity of Variance - raw data
model.levene5 <- leveneTest(model5)
model.levene5  #Result: p-value = 0.6636

#kruskal-wallis test #In my thesis record, I used ANOVA test! I may need to change it to later
model.kruskal5 <-kruskal.test(X17.Dihomo.IsoP ~ day.time, data = metadata)
model.kruskal5  #Kruskal-Wallis chi-squared = 1.4561, df = 3, p-value = 0.6924
dunnTest(X17.Dihomo.IsoP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#6. X10.F4.NeuroP ###############################################################################################################################
#Calculate statistics By day.time
model6 <- lm( X10.F4.NeuroP ~ day.time, metadata) #create a model
anova(model6) #Result: F-value(3,14): 8.4258 ; p-value:  0.001904

#Normality test -
par(mfrow=c(2,2))
plot(model6) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals6 <- residuals(object = model6)
shapiro.test(x = model.residuals6) #result: W = 0.95092, p-value = 0.4396

#Homogeneity of Variance - raw data
model.levene6 <- leveneTest(model6, center = mean)
model.levene6  #Result: p-value = 0.04766

#kruskal-wallis test 
model.kruskal6 <-kruskal.test(X10.F4.NeuroP ~ day.time, data = metadata)
model.kruskal6  #Kruskal-Wallis chi-squared = 11.97, df = 3, p-value = 0.007486
dunnTest(X10.F4.NeuroP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups


#7. X13.F4.NeuroP ###############################################################################################################################
#Calculate statistics By day.time
model7 <- lm( X13.F4.NeuroP ~ day.time, metadata) #create a model
anova(model7) #Result: F-value(3,14): 6.3928 ; p-value:  0.00594 

#Normality test -
par(mfrow=c(2,2))
plot(model7) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals7 <- residuals(object = model7)
shapiro.test(x = model.residuals7) #result: W = 0.9565, p-value = 0.5357

#Homogeneity of Variance - raw data
model.levene7 <- leveneTest(model7, center = mean)
model.levene7  #Result: p-value = 0.005887

#kruskal-wallis test 
model.kruskal7 <-kruskal.test(X13.F4.NeuroP ~ day.time, data = metadata)
model.kruskal7  #Kruskal-Wallis chi-squared = 10.447, df = 3, p-value = 0.01512
dunnTest(X13.F4.NeuroP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#8. X20.F4.NeuroP ###############################################################################################################################
#Calculate statistics By day.time
model8 <- lm( X20.F4.NeuroP ~ day.time, metadata) #create a model
anova(model8) #Result: F-value(3,14): 1.2292 ; p-value:  0.336 

#Normality test -
par(mfrow=c(2,2))
plot(model8) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals8 <- residuals(object = model8)
shapiro.test(x = model.residuals8) #result: W = 0.9565, p-value = 0.5357

#Homogeneity of Variance - raw data
model.levene8 <- leveneTest(model8, center = mean)
model.levene8  #Result: p-value = 0.6318

#9. X9.D1.PhytoP ###############################################################################################################################
#Calculate statistics By day.time
model9 <- lm( X9.D1.PhytoP ~ day.time, metadata) #create a model
anova(model9) #Result: F-value(3,14): 0.5799 ; p-value:  0.6378

#Normality test -
par(mfrow=c(2,2))
plot(model9) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals9 <- residuals(object = model9)
shapiro.test(x = model.residuals9) #result: W = 0.93796, p-value = 0.2673

#Homogeneity of Variance - raw data
model.levene9 <- leveneTest(model9, center = mean)
model.levene9  #Result: p-value = 0.3475

#10. X9.F1.PhytoP ###############################################################################################################################
#Calculate statistics By day.time
model10 <- lm( X9.F1.PhytoP ~ day.time, metadata) #create a model
anova(model10) #Result: F-value(3,14): 1.8873 ; p-value: 0.1781

#Normality test -
par(mfrow=c(2,2))
plot(model10) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals10 <- residuals(object = model10)
shapiro.test(x = model.residuals10) #result: W = 0.92931, p-value = 0.1888

#Homogeneity of Variance - raw data
model.levene10 <- leveneTest(model10, center = mean)
model.levene10  #Result: p-value = 0.000384

#kruskal-wallis test 
model.kruskal10 <-kruskal.test( X9.F1.PhytoP ~ day.time, data = metadata)
model.kruskal10  #Kruskal-Wallis chi-squared = 2.0684, df = 3, p-value = 0.5583
dunnTest( X9.F1.PhytoP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#11. X16.F1.PhytoP ###############################################################################################################################
#Calculate statistics By day.time
model11 <- lm( X16.F1.PhytoP ~ day.time, metadata) #create a model
anova(model11) #Result: F-value(3,14): 3.7217 ; p-value: 0.03706

#Normality test -
par(mfrow=c(2,2))
plot(model11) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals11 <- residuals(object = model11)
shapiro.test(x = model.residuals11) #result: W = 0.90923, p-value = 0.0834

#Homogeneity of Variance - raw data
model.levene11 <- leveneTest(model11, center = mean)
model.levene11  #Result: 0.07278

#Tukey Test
TukeyHSD(aov(X16.F1.PhytoP ~ day.time,metadata))

#12. X16.B1.PhytoP ###############################################################################################################################
#Calculate statistics By day.time
model12 <- lm( X16.B1.PhytoP ~ day.time, metadata) #create a model
anova(model12) #Result: F-value(3,14): 1.5722 ; p-value: 0.2404

#Normality test -
par(mfrow=c(2,2))
plot(model12) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals12 <- residuals(object = model12)
shapiro.test(x = model.residuals12) #result: W = 0.94374, p-value = 0.3353

#Homogeneity of Variance - raw data
model.levene12 <- leveneTest(model12, center = mean)
model.levene12  #Result: p-value = 0.01137

#kruskal-wallis test 
model.kruskal10 <-kruskal.test( X16.B1.PhytoP ~ day.time, data = metadata)
model.kruskal10  #Kruskal-Wallis chi-squared = 3, df = 3, p-value = 0.3916
dunnTest( X16.B1.PhytoP ~ day.time, data = metadata, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#13. X5.F3t.IsoP ###############################################################################################################################
#Calculate statistics By day.time
model13 <- lm( X5.F3t.IsoP ~ day.time, metadata) #create a model
anova(model13) #Result: F-value(3,14): 2.489 ; p-value: 0.103

#Normality test -
par(mfrow=c(2,2))
plot(model13) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals13 <- residuals(object = model13)
shapiro.test(x = model.residuals13) #result: W = 0.89802, p-value = 0.05307

#Homogeneity of Variance - raw data
model.levene13 <- leveneTest(model13, center = mean)
model.levene13  #Result: p-value = 0.3101


#14. X15.F3t.IsoP ###############################################################################################################################
#Calculate statistics By day.time
model14 <- lm( X15.F3t.IsoP ~ day.time, metadata) #create a model
anova(model14) #Result: F-value(3,14): 0.5337 ; p-value:0.6666

#Normality test -
par(mfrow=c(2,2))
plot(model14) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals14 <- residuals(object = model14)
shapiro.test(x = model.residuals14) #result: W = 0.75383, p-value = 0.0003629

#Homogeneity of Variance - raw data
model.levene14 <- leveneTest(model14, center = mean)
model.levene14  #Result: p-value = 0.04671

#kruskal-wallis test 
model.kruskal10 <-kruskal.test( X15.F3t.IsoP ~ day.time, data = metadata)
model.kruskal10  #Kruskal-Wallis chi-squared = 0.56842, df = 3, p-value = 0.9036

#15. X8.F3.IsoP ###############################################################################################################################
#Calculate statistics By day.time
model15 <- lm( X8.F3.IsoP ~ day.time, metadata) #create a model
anova(model15) #Result: F-value(3,15): 0.8128  ; p-value: 0.5078

#Normality test -
par(mfrow=c(2,2))
plot(model15) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

model.residuals15 <- residuals(object = model15)
shapiro.test(x = model.residuals15) #result: W = 0.96915, p-value = 0.7816

#Homogeneity of Variance - raw data
model.levene15 <- leveneTest(model15, center = mean)
model.levene15  #Result: p-value = 0.1305

