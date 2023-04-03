
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity and the species in the genus Allistipes in relation to body weight change
# 2) Study design:  Cross-sectional design
# 3) Endpoints: primary: Long-term body weight change from age 21 to the time of stool collection
#               secondary: BMI and fat mass% at stool collection, body weight change between the two collections of stool sample, and biomarkers of HbA1c and CRP
# 4) Exposures: primary: cumulative average of long-term physical activity level
#               secondary: short-term physical activity level measured by accelerometer, long- and short-term physical activity by intensity
# 5) Covariates: age, diet quality, total energy intake, smoking status, antibiotic use, probiotic use, stool type
# 6) Follow-up: 2012 (Men's Lifestyle Validation Study)  
#############################################################################################################################################

rm(list=ls())
getwd()
setwd("/udd/nhkwa/mlvs")
options(max.print=1000000)

library(viridis)
library(ggplot2)
library(grid)
library(tidyverse)
library(cowplot)
library(data.table)
library(readr)
library(haven)
library(plyr)
library(dplyr)
library(vegan)
library(scales)
library(RColorBrewer)
library(grid)
library(pheatmap)
library(lmerTest)
library(lme4)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(car)
library(Maaslin2)
library(gridExtra)
library(knitr)
library(readxl)
library(ggpubr)
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)
library(Hmisc)

# 8 alistipes species and 1 alistipes genus
#1 s__Alistipes_putredinis
#2 s__Alistipes_finegoldii
#3 s__Alistipes_indistinctus
#4 s__Alistipes_onderdonkii
#5 s__Alistipes_senegalensis
#6 s__Alistipes_shahii
#7 s__Alistipes_sp_AP11
#8 s__Alistipes_sp_HGB5

meta_species<-read.csv(file="./data_generated/meta_species.csv", header = TRUE)
meta_species$mlpa_paqlong <- meta_species$mpa_paqlong+meta_species$lpa_paqlong

myvars <- c("s__Alistipes_putredinis", "s__Alistipes_onderdonkii", "s__Alistipes_shahii", "s__Alistipes_finegoldii", "s__Alistipes_sp_HGB5", "s__Alistipes_indistinctus",  "s__Alistipes_sp_AP11","s__Alistipes_senegalensis")
newdata <- meta_species[myvars]
res2 <- cor(as.matrix(newdata))
res2
res2.cor <- data.frame(res2)
write.csv(res2.cor, file="./data_generated/Allistipes species correlation.csv")

############################
#1. s__Alistipes_putredinis
############################
summary(meta_species$s__Alistipes_putredinis)
table(meta_species$s__Alistipes_putredinis)
graphics.off()
histogram(meta_species$s__Alistipes_putredinis)
meta_species$s__Alistipes_putredinis2c <- with(meta_species, ifelse(s__Alistipes_putredinis> 3.0033, 1, 0))
table(meta_species$s__Alistipes_putredinis2c)
summary(meta_species[meta_species$s__Alistipes_putredinis2c==1,]$s__Alistipes_putredinis)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_putredinis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_putredinis.uni)
bmi_dlw.s__Alistipes_putredinis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_putredinis.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_putredinis.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.below)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.below)
lmer.bmi_dlw.s__Alistipes_putredinis.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.above)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_putredinis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_putredinis.uni)
pfat_dlw.s__Alistipes_putredinis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_putredinis.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_putredinis.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.below)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.below)
lmer.pfat_dlw.s__Alistipes_putredinis.above <- lmer ( pfat_dlw ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.above)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_putredinis.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_putredinis2c + act_paqlong*s__Alistipes_putredinis2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_putredinis.uni)
wtchgsto21.s__Alistipes_putredinis.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_putredinis2c + act_paqlong*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_putredinis.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_putredinis.below <- lmer ( wtchgsto21 ~ act_paqlong
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_putredinis.below)
confint(lmer.wtchgsto21.s__Alistipes_putredinis.below)
lmer.wtchgsto21.s__Alistipes_putredinis.above <- lmer ( wtchgsto21 ~ act_paqlong
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_putredinis.above)
confint(lmer.wtchgsto21.s__Alistipes_putredinis.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_putredinis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_putredinis.uni)
lg2hba1c.s__Alistipes_putredinis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_putredinis.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_putredinis.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.below)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.below)
lmer.lg2hba1c.s__Alistipes_putredinis.above <- lmer ( lg2hba1c ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.above)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.above)

####### lg2crp ########
lg2crp.s__Alistipes_putredinis.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_putredinis.uni)
lg2crp.s__Alistipes_putredinis.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_putredinis.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_putredinis.below <- lmer ( lg2crp ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.below)
confint(lmer.lg2crp.s__Alistipes_putredinis.below)
lmer.lg2crp.s__Alistipes_putredinis.above <- lmer ( lg2crp ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.above)
confint(lmer.lg2crp.s__Alistipes_putredinis.above)

############################
#2. s__Alistipes_finegoldii
############################
summary(meta_species$s__Alistipes_finegoldii)
table(meta_species$s__Alistipes_finegoldii)
graphics.off()
histogram(meta_species$s__Alistipes_finegoldii)
meta_species$s__Alistipes_finegoldii2c <- with(meta_species, ifelse(s__Alistipes_finegoldii> 0.03298, 1, 0))
table(meta_species$s__Alistipes_finegoldii2c)
summary(meta_species[meta_species$s__Alistipes_finegoldii2c==1,]$s__Alistipes_finegoldii)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_finegoldii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_finegoldii.uni)
bmi_dlw.s__Alistipes_finegoldii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_finegoldii.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_finegoldii.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_finegoldii.below)
confint(lmer.bmi_dlw.s__Alistipes_finegoldii.below)
lmer.bmi_dlw.s__Alistipes_finegoldii.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_finegoldii.above)
confint(lmer.bmi_dlw.s__Alistipes_finegoldii.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_finegoldii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_finegoldii.uni)
pfat_dlw.s__Alistipes_finegoldii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_finegoldii.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_finegoldii.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_finegoldii.below)
confint(lmer.pfat_dlw.s__Alistipes_finegoldii.below)
lmer.pfat_dlw.s__Alistipes_finegoldii.above <- lmer ( pfat_dlw ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_finegoldii.above)
confint(lmer.pfat_dlw.s__Alistipes_finegoldii.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_finegoldii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_finegoldii2c + act_paqlong*s__Alistipes_finegoldii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_finegoldii.uni)
wtchgsto21.s__Alistipes_finegoldii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_finegoldii2c + act_paqlong*s__Alistipes_finegoldii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_finegoldii.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_finegoldii.below <- lmer ( wtchgsto21 ~ act_paqlong 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_finegoldii.below)
confint(lmer.wtchgsto21.s__Alistipes_finegoldii.below)
lmer.wtchgsto21.s__Alistipes_finegoldii.above <- lmer ( wtchgsto21 ~ act_paqlong
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_finegoldii.above)
confint(lmer.wtchgsto21.s__Alistipes_finegoldii.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_finegoldii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_finegoldii.uni)
lg2hba1c.s__Alistipes_finegoldii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_finegoldii.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_finegoldii.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_finegoldii.below)
confint(lmer.lg2hba1c.s__Alistipes_finegoldii.below)
lmer.lg2hba1c.s__Alistipes_finegoldii.above <- lmer ( lg2hba1c ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_finegoldii.above)
confint(lmer.lg2hba1c.s__Alistipes_finegoldii.above)

####### lg2crp ########
lg2crp.s__Alistipes_finegoldii.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_finegoldii.uni)
lg2crp.s__Alistipes_finegoldii.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_finegoldii.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_finegoldii.below <- lmer ( lg2crp ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==0, ])
summary(lmer.lg2crp.s__Alistipes_finegoldii.below)
confint(lmer.lg2crp.s__Alistipes_finegoldii.below)
lmer.lg2crp.s__Alistipes_finegoldii.above <- lmer ( lg2crp ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_finegoldii2c==1, ])
summary(lmer.lg2crp.s__Alistipes_finegoldii.above)
confint(lmer.lg2crp.s__Alistipes_finegoldii.above)

##############################
#3. s__Alistipes_indistinctus
##############################
summary(meta_species$s__Alistipes_indistinctus)
table(meta_species$s__Alistipes_indistinctus)
graphics.off()
histogram(meta_species$s__Alistipes_indistinctus)
meta_species$s__Alistipes_indistinctus2c <- with(meta_species, ifelse(s__Alistipes_indistinctus> 0, 1, 0))
table(meta_species$s__Alistipes_indistinctus2c)
summary(meta_species[meta_species$s__Alistipes_indistinctus2c==1,]$s__Alistipes_indistinctus)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_indistinctus.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_indistinctus.uni)
bmi_dlw.s__Alistipes_indistinctus.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_indistinctus.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_indistinctus.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_indistinctus.below)
confint(lmer.bmi_dlw.s__Alistipes_indistinctus.below)
lmer.bmi_dlw.s__Alistipes_indistinctus.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_indistinctus.above)
confint(lmer.bmi_dlw.s__Alistipes_indistinctus.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_indistinctus.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_indistinctus.uni)
pfat_dlw.s__Alistipes_indistinctus.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_indistinctus.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_indistinctus.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_indistinctus.below)
confint(lmer.pfat_dlw.s__Alistipes_indistinctus.below)
lmer.pfat_dlw.s__Alistipes_indistinctus.above <- lmer ( pfat_dlw ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_indistinctus.above)
confint(lmer.pfat_dlw.s__Alistipes_indistinctus.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_indistinctus.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_indistinctus2c + act_paqlong*s__Alistipes_indistinctus2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_indistinctus.uni)
wtchgsto21.s__Alistipes_indistinctus.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_indistinctus2c + act_paqlong*s__Alistipes_indistinctus2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_indistinctus.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_indistinctus.below <- lmer ( wtchgsto21 ~ act_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_indistinctus.below)
confint(lmer.wtchgsto21.s__Alistipes_indistinctus.below)
lmer.wtchgsto21.s__Alistipes_indistinctus.above <- lmer ( wtchgsto21 ~ act_paqlong
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_indistinctus.above)
confint(lmer.wtchgsto21.s__Alistipes_indistinctus.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_indistinctus.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_indistinctus.uni)
lg2hba1c.s__Alistipes_indistinctus.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_indistinctus.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_indistinctus.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_indistinctus.below)
confint(lmer.lg2hba1c.s__Alistipes_indistinctus.below)
lmer.lg2hba1c.s__Alistipes_indistinctus.above <- lmer ( lg2hba1c ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_indistinctus.above)
confint(lmer.lg2hba1c.s__Alistipes_indistinctus.above)

####### lg2crp ########
lg2crp.s__Alistipes_indistinctus.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_indistinctus.uni)
lg2crp.s__Alistipes_indistinctus.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_indistinctus.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_indistinctus.below <- lmer ( lg2crp ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==0, ])
summary(lmer.lg2crp.s__Alistipes_indistinctus.below)
confint(lmer.lg2crp.s__Alistipes_indistinctus.below)
lmer.lg2crp.s__Alistipes_indistinctus.above <- lmer ( lg2crp ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_indistinctus2c==1, ])
summary(lmer.lg2crp.s__Alistipes_indistinctus.above)
confint(lmer.lg2crp.s__Alistipes_indistinctus.above)

#############################
#4. s__Alistipes_onderdonkii
#############################
summary(meta_species$s__Alistipes_onderdonkii)
table(meta_species$s__Alistipes_onderdonkii)
graphics.off()
histogram(meta_species$s__Alistipes_onderdonkii)
meta_species$s__Alistipes_onderdonkii2c <- with(meta_species, ifelse(s__Alistipes_onderdonkii> 0.43955, 1, 0))
table(meta_species$s__Alistipes_onderdonkii2c)
summary(meta_species[meta_species$s__Alistipes_onderdonkii2c==1,]$s__Alistipes_onderdonkii)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_onderdonkii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_onderdonkii.uni)
bmi_dlw.s__Alistipes_onderdonkii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_onderdonkii.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_onderdonkii.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_onderdonkii.below)
confint(lmer.bmi_dlw.s__Alistipes_onderdonkii.below)
lmer.bmi_dlw.s__Alistipes_onderdonkii.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_onderdonkii.above)
confint(lmer.bmi_dlw.s__Alistipes_onderdonkii.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_onderdonkii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_onderdonkii.uni)
pfat_dlw.s__Alistipes_onderdonkii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_onderdonkii.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_onderdonkii.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_onderdonkii.below)
confint(lmer.pfat_dlw.s__Alistipes_onderdonkii.below)
lmer.pfat_dlw.s__Alistipes_onderdonkii.above <- lmer ( pfat_dlw ~ paee_pamwk
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_onderdonkii.above)
confint(lmer.pfat_dlw.s__Alistipes_onderdonkii.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_onderdonkii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_onderdonkii2c + act_paqlong*s__Alistipes_onderdonkii2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_onderdonkii.uni)
wtchgsto21.s__Alistipes_onderdonkii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_onderdonkii2c + act_paqlong*s__Alistipes_onderdonkii2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_onderdonkii.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_onderdonkii.below <- lmer ( wtchgsto21 ~ act_paqlong 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_onderdonkii.below)
confint(lmer.wtchgsto21.s__Alistipes_onderdonkii.below)
lmer.wtchgsto21.s__Alistipes_onderdonkii.above <- lmer ( wtchgsto21 ~ act_paqlong
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_onderdonkii.above)
confint(lmer.wtchgsto21.s__Alistipes_onderdonkii.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_onderdonkii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_onderdonkii.uni)
lg2hba1c.s__Alistipes_onderdonkii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_onderdonkii.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_onderdonkii.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_onderdonkii.below)
confint(lmer.lg2hba1c.s__Alistipes_onderdonkii.below)
lmer.lg2hba1c.s__Alistipes_onderdonkii.above <- lmer ( lg2hba1c ~ paee_pamwk
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_onderdonkii.above)
confint(lmer.lg2hba1c.s__Alistipes_onderdonkii.above)

####### lg2crp ########
lg2crp.s__Alistipes_onderdonkii.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_onderdonkii.uni)
lg2crp.s__Alistipes_onderdonkii.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_onderdonkii.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_onderdonkii.below <- lmer ( lg2crp ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==0, ])
summary(lmer.lg2crp.s__Alistipes_onderdonkii.below)
confint(lmer.lg2crp.s__Alistipes_onderdonkii.below)
lmer.lg2crp.s__Alistipes_onderdonkii.above <- lmer ( lg2crp ~ paee_pamwk
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_onderdonkii2c==1, ])
summary(lmer.lg2crp.s__Alistipes_onderdonkii.above)
confint(lmer.lg2crp.s__Alistipes_onderdonkii.above)

##############################
#5. s__Alistipes_senegalensis
##############################
summary(meta_species$s__Alistipes_senegalensis)
table(meta_species$s__Alistipes_senegalensis)
graphics.off()
histogram(meta_species$s__Alistipes_senegalensis)
meta_species$s__Alistipes_senegalensis2c <- with(meta_species, ifelse(s__Alistipes_senegalensis> 0, 1, 0))
table(meta_species$s__Alistipes_senegalensis2c)
summary(meta_species[meta_species$s__Alistipes_senegalensis2c==1,]$s__Alistipes_senegalensis)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_senegalensis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_senegalensis.uni)
bmi_dlw.s__Alistipes_senegalensis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_senegalensis.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_senegalensis.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_senegalensis.below)
confint(lmer.bmi_dlw.s__Alistipes_senegalensis.below)
lmer.bmi_dlw.s__Alistipes_senegalensis.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_senegalensis.above)
confint(lmer.bmi_dlw.s__Alistipes_senegalensis.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_senegalensis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_senegalensis.uni)
pfat_dlw.s__Alistipes_senegalensis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_senegalensis.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_senegalensis.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_senegalensis.below)
confint(lmer.pfat_dlw.s__Alistipes_senegalensis.below)
lmer.pfat_dlw.s__Alistipes_senegalensis.above <- lmer ( pfat_dlw ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_senegalensis.above)
confint(lmer.pfat_dlw.s__Alistipes_senegalensis.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_senegalensis.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_senegalensis2c + act_paqlong*s__Alistipes_senegalensis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_senegalensis.uni)
wtchgsto21.s__Alistipes_senegalensis.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_senegalensis2c + act_paqlong*s__Alistipes_senegalensis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_senegalensis.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_senegalensis.below <- lmer ( wtchgsto21 ~ act_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_senegalensis.below)
confint(lmer.wtchgsto21.s__Alistipes_senegalensis.below)
lmer.wtchgsto21.s__Alistipes_senegalensis.above <- lmer ( wtchgsto21 ~ act_paqlong
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_senegalensis.above)
confint(lmer.wtchgsto21.s__Alistipes_senegalensis.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_senegalensis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_senegalensis.uni)
lg2hba1c.s__Alistipes_senegalensis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_senegalensis.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_senegalensis.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_senegalensis.below)
confint(lmer.lg2hba1c.s__Alistipes_senegalensis.below)
lmer.lg2hba1c.s__Alistipes_senegalensis.above <- lmer ( lg2hba1c ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_senegalensis.above)
confint(lmer.lg2hba1c.s__Alistipes_senegalensis.above)

####### lg2crp ########
lg2crp.s__Alistipes_senegalensis.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_senegalensis.uni)
lg2crp.s__Alistipes_senegalensis.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_senegalensis.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_senegalensis.below <- lmer ( lg2crp ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==0, ])
summary(lmer.lg2crp.s__Alistipes_senegalensis.below)
confint(lmer.lg2crp.s__Alistipes_senegalensis.below)
lmer.lg2crp.s__Alistipes_senegalensis.above <- lmer ( lg2crp ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_senegalensis2c==1, ])
summary(lmer.lg2crp.s__Alistipes_senegalensis.above)
confint(lmer.lg2crp.s__Alistipes_senegalensis.above)

########################
#6. s__Alistipes_shahii
########################
summary(meta_species$s__Alistipes_shahii)
table(meta_species$s__Alistipes_shahii)
graphics.off()
histogram(meta_species$s__Alistipes_shahii)
meta_species$s__Alistipes_shahii2c <- with(meta_species, ifelse(s__Alistipes_shahii> 0.6692, 1, 0))
table(meta_species$s__Alistipes_shahii2c)
summary(meta_species[meta_species$s__Alistipes_shahii2c==1,]$s__Alistipes_shahii)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_shahii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_shahii.uni)
bmi_dlw.s__Alistipes_shahii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_shahii.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_shahii.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_shahii.below)
confint(lmer.bmi_dlw.s__Alistipes_shahii.below)
lmer.bmi_dlw.s__Alistipes_shahii.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_shahii.above)
confint(lmer.bmi_dlw.s__Alistipes_shahii.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_shahii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                        + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_shahii.uni)
pfat_dlw.s__Alistipes_shahii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                        +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_shahii.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_shahii.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_shahii.below)
confint(lmer.pfat_dlw.s__Alistipes_shahii.below)
lmer.pfat_dlw.s__Alistipes_shahii.above <- lmer ( pfat_dlw ~ paee_pamwk
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_shahii.above)
confint(lmer.pfat_dlw.s__Alistipes_shahii.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_shahii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_shahii2c + act_paqlong*s__Alistipes_shahii2c 
                                        + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_shahii.uni)
wtchgsto21.s__Alistipes_shahii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_shahii2c + act_paqlong*s__Alistipes_shahii2c 
                                        +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_shahii.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_shahii.below <- lmer ( wtchgsto21 ~ act_paqlong 
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_shahii.below)
confint(lmer.wtchgsto21.s__Alistipes_shahii.below)
lmer.wtchgsto21.s__Alistipes_shahii.above <- lmer ( wtchgsto21 ~ act_paqlong
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_shahii.above)
confint(lmer.wtchgsto21.s__Alistipes_shahii.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_shahii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                        + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_shahii.uni)
lg2hba1c.s__Alistipes_shahii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                        +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_shahii.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_shahii.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_shahii.below)
confint(lmer.lg2hba1c.s__Alistipes_shahii.below)
lmer.lg2hba1c.s__Alistipes_shahii.above <- lmer ( lg2hba1c ~ paee_pamwk
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_shahii.above)
confint(lmer.lg2hba1c.s__Alistipes_shahii.above)

####### lg2crp ########
lg2crp.s__Alistipes_shahii.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                        + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_shahii.uni)
lg2crp.s__Alistipes_shahii.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                        +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_shahii.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_shahii.below <- lmer ( lg2crp ~ paee_pamwk 
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==0, ])
summary(lmer.lg2crp.s__Alistipes_shahii.below)
confint(lmer.lg2crp.s__Alistipes_shahii.below)
lmer.lg2crp.s__Alistipes_shahii.above <- lmer ( lg2crp ~ paee_pamwk
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_shahii2c==1, ])
summary(lmer.lg2crp.s__Alistipes_shahii.above)
confint(lmer.lg2crp.s__Alistipes_shahii.above)

#########################
#7. s__Alistipes_sp_AP11
#########################
summary(meta_species$s__Alistipes_sp_AP11)
table(meta_species$s__Alistipes_sp_AP11)
graphics.off()
histogram(meta_species$s__Alistipes_sp_AP11)
meta_species$s__Alistipes_sp_AP112c <- with(meta_species, ifelse(s__Alistipes_sp_AP11> 0, 1, 0))
table(meta_species$s__Alistipes_sp_AP112c)
summary(meta_species[meta_species$s__Alistipes_sp_AP112c==1,]$s__Alistipes_sp_AP11)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_sp_AP11.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                        + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_sp_AP11.uni)

bmi_dlw.s__Alistipes_sp_AP11.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                        +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_sp_AP11.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_sp_AP11.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_sp_AP11.below)
confint(lmer.bmi_dlw.s__Alistipes_sp_AP11.below)
lmer.bmi_dlw.s__Alistipes_sp_AP11.above <- lmer ( bmi_dlw ~ paee_pamwk
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_sp_AP11.above)
confint(lmer.bmi_dlw.s__Alistipes_sp_AP11.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_sp_AP11.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                         + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_sp_AP11.uni)
pfat_dlw.s__Alistipes_sp_AP11.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_sp_AP11.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_sp_AP11.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_sp_AP11.below)
confint(lmer.pfat_dlw.s__Alistipes_sp_AP11.below)
lmer.pfat_dlw.s__Alistipes_sp_AP11.above <- lmer ( pfat_dlw ~ paee_pamwk
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_sp_AP11.above)
confint(lmer.pfat_dlw.s__Alistipes_sp_AP11.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_sp_AP11.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_sp_AP112c + act_paqlong*s__Alistipes_sp_AP112c 
                                         + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_sp_AP11.uni)
wtchgsto21.s__Alistipes_sp_AP11.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_sp_AP112c + act_paqlong*s__Alistipes_sp_AP112c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_sp_AP11.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_sp_AP11.below <- lmer ( wtchgsto21 ~ act_paqlong 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_sp_AP11.below)
confint(lmer.wtchgsto21.s__Alistipes_sp_AP11.below)
lmer.wtchgsto21.s__Alistipes_sp_AP11.above <- lmer ( wtchgsto21 ~ act_paqlong
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_sp_AP11.above)
confint(lmer.wtchgsto21.s__Alistipes_sp_AP11.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_sp_AP11.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                         + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_sp_AP11.uni)
lg2hba1c.s__Alistipes_sp_AP11.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_sp_AP11.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_sp_AP11.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_sp_AP11.below)
confint(lmer.lg2hba1c.s__Alistipes_sp_AP11.below)
lmer.lg2hba1c.s__Alistipes_sp_AP11.above <- lmer ( lg2hba1c ~ paee_pamwk
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_sp_AP11.above)
confint(lmer.lg2hba1c.s__Alistipes_sp_AP11.above)

####### lg2crp ########
lg2crp.s__Alistipes_sp_AP11.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                         + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_sp_AP11.uni)
lg2crp.s__Alistipes_sp_AP11.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_sp_AP11.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_sp_AP11.below <- lmer ( lg2crp ~ paee_pamwk 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==0, ])
summary(lmer.lg2crp.s__Alistipes_sp_AP11.below)
confint(lmer.lg2crp.s__Alistipes_sp_AP11.below)
lmer.lg2crp.s__Alistipes_sp_AP11.above <- lmer ( lg2crp ~ paee_pamwk
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_AP112c==1, ])
summary(lmer.lg2crp.s__Alistipes_sp_AP11.above)
confint(lmer.lg2crp.s__Alistipes_sp_AP11.above)

##########################
#8. s__Alistipes_sp_HGB5
##########################
summary(meta_species$s__Alistipes_sp_HGB5)
table(meta_species$s__Alistipes_sp_HGB5)
graphics.off()
histogram(meta_species$s__Alistipes_sp_HGB5)
meta_species$s__Alistipes_sp_HGB52c <- with(meta_species, ifelse(s__Alistipes_sp_HGB5> 0, 1, 0))
table(meta_species$s__Alistipessp_HGB52c)
summary(meta_species[meta_species$s__Alistipes_sp_HGB52c==1,]$s__Alistipes_sp_HGB5)

####### bmi_dlw ########
bmi_dlw.s__Alistipes_sp_HGB5.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_sp_HGB5.uni)
bmi_dlw.s__Alistipes_sp_HGB5.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_sp_HGB5.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_sp_HGB5.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_sp_HGB5.below)
confint(lmer.bmi_dlw.s__Alistipes_sp_HGB5.below)
lmer.bmi_dlw.s__Alistipes_sp_HGB5.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_sp_HGB5.above)
confint(lmer.bmi_dlw.s__Alistipes_sp_HGB5.above)

####### pfat_dlw ########
pfat_dlw.s__Alistipes_sp_HGB5.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_sp_HGB5.uni)
pfat_dlw.s__Alistipes_sp_HGB5.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_sp_HGB5.uni)

#all act in each stratum
lmer.pfat_dlw.s__Alistipes_sp_HGB5.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_sp_HGB5.below)
confint(lmer.pfat_dlw.s__Alistipes_sp_HGB5.below)
lmer.pfat_dlw.s__Alistipes_sp_HGB5.above <- lmer ( pfat_dlw ~ paee_pamwk
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_sp_HGB5.above)
confint(lmer.pfat_dlw.s__Alistipes_sp_HGB5.above)

####### wtchgsto21 ########
wtchgsto21.s__Alistipes_sp_HGB5.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_sp_HGB52c + act_paqlong*s__Alistipes_sp_HGB52c 
                                         + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_sp_HGB5.uni)
wtchgsto21.s__Alistipes_sp_HGB5.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Alistipes_sp_HGB52c + act_paqlong*s__Alistipes_sp_HGB52c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.s__Alistipes_sp_HGB5.uni)

#all act in each stratum
lmer.wtchgsto21.s__Alistipes_sp_HGB5.below <- lmer ( wtchgsto21 ~ act_paqlong 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==0, ])
summary(lmer.wtchgsto21.s__Alistipes_sp_HGB5.below)
confint(lmer.wtchgsto21.s__Alistipes_sp_HGB5.below)
lmer.wtchgsto21.s__Alistipes_sp_HGB5.above <- lmer ( wtchgsto21 ~ act_paqlong
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==1, ])
summary(lmer.wtchgsto21.s__Alistipes_sp_HGB5.above)
confint(lmer.wtchgsto21.s__Alistipes_sp_HGB5.above)

####### lg2hba1c ########
lg2hba1c.s__Alistipes_sp_HGB5.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_sp_HGB5.uni)
lg2hba1c.s__Alistipes_sp_HGB5.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_sp_HGB5.uni)

#all act in each stratum
lmer.lg2hba1c.s__Alistipes_sp_HGB5.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_sp_HGB5.below)
confint(lmer.lg2hba1c.s__Alistipes_sp_HGB5.below)
lmer.lg2hba1c.s__Alistipes_sp_HGB5.above <- lmer ( lg2hba1c ~ paee_pamwk
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_sp_HGB5.above)
confint(lmer.lg2hba1c.s__Alistipes_sp_HGB5.above)

####### lg2crp ########
lg2crp.s__Alistipes_sp_HGB5.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_sp_HGB5.uni)
lg2crp.s__Alistipes_sp_HGB5.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_sp_HGB5.uni)

#all act in each stratum
lmer.lg2crp.s__Alistipes_sp_HGB5.below <- lmer ( lg2crp ~ paee_pamwk 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==0, ])
summary(lmer.lg2crp.s__Alistipes_sp_HGB5.below)
confint(lmer.lg2crp.s__Alistipes_sp_HGB5.below)
lmer.lg2crp.s__Alistipes_sp_HGB5.above <- lmer ( lg2crp ~ paee_pamwk
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_sp_HGB52c==1, ])
summary(lmer.lg2crp.s__Alistipes_sp_HGB5.above)
confint(lmer.lg2crp.s__Alistipes_sp_HGB5.above)
