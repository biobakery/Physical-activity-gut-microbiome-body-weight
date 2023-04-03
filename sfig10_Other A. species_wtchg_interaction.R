
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity measured at the 1st stool collection and the species in the genus Allistipes in relation to body weight change between the. 2 stool samples
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

meta_species<-read.csv(file="./data_generated/meta_species.csv", header = TRUE)
meta_species_12stool <- subset(meta_species, cvisit==1)

####### weightchg_blood ########
############################
#1. s__Alistipes_putredinis
############################
summary(meta_species_12stool$s__Alistipes_putredinis)
table(meta_species_12stool$s__Alistipes_putredinis)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_putredinis)
meta_species_12stool$s__Alistipes_putredinis2c <- with(meta_species_12stool, ifelse(s__Alistipes_putredinis>2.9085 , 1, 0))
table(meta_species_12stool$s__Alistipes_putredinis2c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1,]$s__Alistipes_putredinis)

weightchg_blood.s__Alistipes_putredinis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_putredinis.uni)
weightchg_blood.s__Alistipes_putredinis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_putredinis.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_putredinis.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.below)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.below)
lmer.weightchg_blood.s__Alistipes_putredinis.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.above)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.above)

############################
#2. s__Alistipes_finegoldii
############################
summary(meta_species_12stool$s__Alistipes_finegoldii)
table(meta_species_12stool$s__Alistipes_finegoldii)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_finegoldii)
meta_species_12stool$s__Alistipes_finegoldii2c <- with(meta_species_12stool, ifelse(s__Alistipes_finegoldii>0.03458 , 1, 0))
table(meta_species_12stool$s__Alistipes_finegoldii2c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_finegoldii2c==1,]$s__Alistipes_finegoldii)
weightchg_blood.s__Alistipes_finegoldii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_finegoldii.uni)
weightchg_blood.s__Alistipes_finegoldii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_finegoldii2c + paee_pamwk*s__Alistipes_finegoldii2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_finegoldii.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_finegoldii.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_finegoldii2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_finegoldii.below)
confint(lmer.weightchg_blood.s__Alistipes_finegoldii.below)
lmer.weightchg_blood.s__Alistipes_finegoldii.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_finegoldii2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_finegoldii.above)
confint(lmer.weightchg_blood.s__Alistipes_finegoldii.above)

##############################
#3. s__Alistipes_indistinctus
##############################
summary(meta_species_12stool$s__Alistipes_indistinctus)
table(meta_species_12stool$s__Alistipes_indistinctus)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_indistinctus)
meta_species_12stool$s__Alistipes_indistinctus2c <- with(meta_species_12stool, ifelse(s__Alistipes_indistinctus>0.001735 , 1, 0))
table(meta_species_12stool$s__Alistipes_indistinctus2c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_indistinctus2c==1,]$s__Alistipes_indistinctus)
weightchg_blood.s__Alistipes_indistinctus.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_indistinctus.uni)
weightchg_blood.s__Alistipes_indistinctus.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_indistinctus2c + paee_pamwk*s__Alistipes_indistinctus2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_indistinctus.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_indistinctus.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_indistinctus2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_indistinctus.below)
confint(lmer.weightchg_blood.s__Alistipes_indistinctus.below)
lmer.weightchg_blood.s__Alistipes_indistinctus.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_indistinctus2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_indistinctus.above)
confint(lmer.weightchg_blood.s__Alistipes_indistinctus.above)

#############################
#4. s__Alistipes_onderdonkii
#############################
summary(meta_species_12stool$s__Alistipes_onderdonkii)
table(meta_species_12stool$s__Alistipes_onderdonkii)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_onderdonkii)
meta_species_12stool$s__Alistipes_onderdonkii2c <- with(meta_species_12stool, ifelse(s__Alistipes_onderdonkii>0.42579 , 1, 0))
table(meta_species_12stool$s__Alistipes_onderdonkii2c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_onderdonkii2c==1,]$s__Alistipes_onderdonkii)
weightchg_blood.s__Alistipes_onderdonkii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_onderdonkii.uni)
weightchg_blood.s__Alistipes_onderdonkii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_onderdonkii2c + paee_pamwk*s__Alistipes_onderdonkii2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_onderdonkii.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_onderdonkii.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_onderdonkii2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_onderdonkii.below)
confint(lmer.weightchg_blood.s__Alistipes_onderdonkii.below)
lmer.weightchg_blood.s__Alistipes_onderdonkii.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_onderdonkii2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_onderdonkii.above)
confint(lmer.weightchg_blood.s__Alistipes_onderdonkii.above)

##############################
#5. s__Alistipes_senegalensis
##############################
summary(meta_species_12stool$s__Alistipes_senegalensis)
table(meta_species_12stool$s__Alistipes_senegalensis)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_senegalensis)
meta_species_12stool$s__Alistipes_senegalensis2c <- with(meta_species_12stool, ifelse(s__Alistipes_senegalensis>0 , 1, 0))
table(meta_species_12stool$s__Alistipes_senegalensis2c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_senegalensis2c==1,]$s__Alistipes_senegalensis)
weightchg_blood.s__Alistipes_senegalensis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_senegalensis.uni)
weightchg_blood.s__Alistipes_senegalensis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_senegalensis2c + paee_pamwk*s__Alistipes_senegalensis2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_senegalensis.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_senegalensis.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_senegalensis2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_senegalensis.below)
confint(lmer.weightchg_blood.s__Alistipes_senegalensis.below)
lmer.weightchg_blood.s__Alistipes_senegalensis.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_senegalensis2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_senegalensis.above)
confint(lmer.weightchg_blood.s__Alistipes_senegalensis.above)

########################
#6. s__Alistipes_shahii
########################
summary(meta_species_12stool$s__Alistipes_shahii)
table(meta_species_12stool$s__Alistipes_shahii)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_shahii)
meta_species_12stool$s__Alistipes_shahii2c <- with(meta_species_12stool, ifelse(s__Alistipes_shahii>0.7348 , 1, 0))
table(meta_species_12stool$s__Alistipes_shahii2c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_shahii2c==1,]$s__Alistipes_shahii)
weightchg_blood.s__Alistipes_shahii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_shahii.uni)
weightchg_blood.s__Alistipes_shahii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_shahii2c + paee_pamwk*s__Alistipes_shahii2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_shahii.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_shahii.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_shahii2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_shahii.below)
confint(lmer.weightchg_blood.s__Alistipes_shahii.below)
lmer.weightchg_blood.s__Alistipes_shahii.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_shahii2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_shahii.above)
confint(lmer.weightchg_blood.s__Alistipes_shahii.above)

#########################
#7. s__Alistipes_sp_AP11
#########################
summary(meta_species_12stool$s__Alistipes_sp_AP11)
table(meta_species_12stool$s__Alistipes_sp_AP11)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_sp_AP11)
meta_species_12stool$s__Alistipes_sp_AP112c <- with(meta_species_12stool, ifelse(s__Alistipes_sp_AP11>0 , 1, 0))
table(meta_species_12stool$s__Alistipes_sp_AP112c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_sp_AP112c==1,]$s__Alistipes_sp_AP11)
weightchg_blood.s__Alistipes_sp_AP11.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_sp_AP11.uni)
weightchg_blood.s__Alistipes_sp_AP11.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_sp_AP112c + paee_pamwk*s__Alistipes_sp_AP112c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_sp_AP11.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_sp_AP11.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_sp_AP112c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_sp_AP11.below)
confint(lmer.weightchg_blood.s__Alistipes_sp_AP11.below)
lmer.weightchg_blood.s__Alistipes_sp_AP11.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_sp_AP112c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_sp_AP11.above)
confint(lmer.weightchg_blood.s__Alistipes_sp_AP11.above)

##########################
#8. s__Alistipes_sp_HGB5
##########################
summary(meta_species_12stool$s__Alistipes_sp_HGB5)
table(meta_species_12stool$s__Alistipes_sp_HGB5)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_sp_HGB5)
meta_species_12stool$s__Alistipes_sp_HGB52c <- with(meta_species_12stool, ifelse(s__Alistipes_sp_HGB5>0 , 1, 0))
table(meta_species_12stool$s__Alistipes_sp_HGB52c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_sp_HGB52c==1,]$s__Alistipes_sp_HGB5)

weightchg_blood.s__Alistipes_sp_HGB5.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_sp_HGB5.uni)
weightchg_blood.s__Alistipes_sp_HGB5.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_sp_HGB52c + paee_pamwk*s__Alistipes_sp_HGB52c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_sp_HGB5.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_sp_HGB5.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_sp_HGB52c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_sp_HGB5.below)
confint(lmer.weightchg_blood.s__Alistipes_sp_HGB5.below)
lmer.weightchg_blood.s__Alistipes_sp_HGB5.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_sp_HGB52c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_sp_HGB5.above)
confint(lmer.weightchg_blood.s__Alistipes_sp_HGB5.above)
