
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity measured at the 1st stool collection and the top 10 abundant species in relation to body weight change between the 2 collections of stoll samples
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
meta_species$paee_pamwk <- meta_species$paee_pam*7
meta_species_12stool <- subset(meta_species, cvisit==1)

# top 10 most abundant species
# 1. s__Eubacterium_rectale                             
# 2. s__Faecalibacterium_prausnitzii
# 3. s__Subdoligranulum_unclassified
# 4. s__Bacteroides_uniformis                                                 
# 5. s__Prevotella_copri
# 6. s__Alistipes_putredinis
# 7. s__Eubacterium_siraeum                                              
# 8. s__Ruminococcus_bromii
# 9. s__Bacteroides_vulgatus
# 10.s__Bacteroides_stercoris

####### weightchg_blood and % ########

#1. s__Eubacterium_rectale
summary(meta_species_12stool$s__Eubacterium_rectale)
table(meta_species_12stool$s__Eubacterium_rectale)
graphics.off()
histogram(meta_species_12stool$s__Eubacterium_rectale)
meta_species_12stool$s__Eubacterium_rectale2c <- with(meta_species_12stool, ifelse(s__Eubacterium_rectale> 6.518, 1, 0))
table(meta_species_12stool$s__Eubacterium_rectale2c)
summary(meta_species_12stool[meta_species_12stool$s__Eubacterium_rectale2c==1,]$s__Eubacterium_rectale)
#chagne
weightchg_blood.s__Eubacterium_rectale.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                          + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Eubacterium_rectale.uni)

weightchg_blood.s__Eubacterium_rectale.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                                  +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Eubacterium_rectale.uni)
#chagne %
weightchgp_blood.s__Eubacterium_rectale.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                                   + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Eubacterium_rectale.uni)

weightchgp_blood.s__Eubacterium_rectale.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Eubacterium_rectale.uni)

#2. s__Faecalibacterium_prausnitzii
summary(meta_species_12stool$s__Faecalibacterium_prausnitzii)
table(meta_species_12stool$s__Faecalibacterium_prausnitzii)
graphics.off()
histogram(meta_species_12stool$s__Faecalibacterium_prausnitzii)
meta_species_12stool$s__Faecalibacterium_prausnitzii2c <- with(meta_species_12stool, ifelse(s__Faecalibacterium_prausnitzii> 6.452, 1, 0))
table(meta_species_12stool$s__Faecalibacterium_prausnitzii2c)
summary(meta_species_12stool[meta_species_12stool$s__Faecalibacterium_prausnitzii2c==1,]$s__Faecalibacterium_prausnitzii)
#chagne
weightchg_blood.s__Faecalibacterium_prausnitzii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                  + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Faecalibacterium_prausnitzii.uni)

weightchg_blood.s__Faecalibacterium_prausnitzii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Faecalibacterium_prausnitzii.uni)
#chagne %
weightchgp_blood.s__Faecalibacterium_prausnitzii.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                            + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Faecalibacterium_prausnitzii.uni)

weightchgp_blood.s__Faecalibacterium_prausnitzii.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Faecalibacterium_prausnitzii.uni)

#3. s__Subdoligranulum_unclassified
summary(meta_species_12stool$s__Subdoligranulum_unclassified)
table(meta_species_12stool$s__Subdoligranulum_unclassified)
graphics.off()
histogram(meta_species_12stool$s__Subdoligranulum_unclassified)
meta_species_12stool$s__Subdoligranulum_unclassified2c <- with(meta_species_12stool, ifelse(s__Subdoligranulum_unclassified> 5.037, 1, 0))
table(meta_species_12stool$s__Subdoligranulum_unclassified2c)
summary(meta_species_12stool[meta_species_12stool$s__Subdoligranulum_unclassified2c==1,]$s__Subdoligranulum_unclassified)
#chagne

weightchg_blood.s__Subdoligranulum_unclassified.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Subdoligranulum_unclassified.uni)

weightchg_blood.s__Subdoligranulum_unclassified.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Subdoligranulum_unclassified.uni)
#chagne %
weightchgp_blood.s__Subdoligranulum_unclassified.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                            + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Subdoligranulum_unclassified.uni)

weightchgp_blood.s__Subdoligranulum_unclassified.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Subdoligranulum_unclassified.uni)

#4. s__Bacteroides_uniformis
summary(meta_species_12stool$s__Bacteroides_uniformis)
table(meta_species_12stool$s__Bacteroides_uniformis)
graphics.off()
histogram(meta_species_12stool$s__Bacteroides_uniformis)
meta_species_12stool$s__Bacteroides_uniformis2c <- with(meta_species_12stool, ifelse(s__Bacteroides_uniformis> 4.059, 1, 0))
table(meta_species_12stool$s__Bacteroides_uniformis2c)
summary(meta_species_12stool[meta_species_12stool$s__Bacteroides_uniformis2c==1,]$s__Bacteroides_uniformis)
#chagne
weightchg_blood.s__Bacteroides_uniformis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Bacteroides_uniformis.uni)

weightchg_blood.s__Bacteroides_uniformis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Bacteroides_uniformis.uni)
#chagne %
weightchgp_blood.s__Bacteroides_uniformis.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                                     + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Bacteroides_uniformis.uni)

weightchgp_blood.s__Bacteroides_uniformis.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Bacteroides_uniformis.uni)

#5. s__Prevotella_copri
summary(meta_species_12stool$s__Prevotella_copri)
table(meta_species_12stool$s__Prevotella_copri)
graphics.off()
histogram(meta_species_12stool$s__Prevotella_copri)
meta_species_12stool$s__Prevotella_copri2c <- with(meta_species_12stool, ifelse(s__Prevotella_copri> 0, 1, 0))
table(meta_species_12stool$s__Prevotella_copri2c)
summary(meta_species_12stool[meta_species_12stool$s__Prevotella_copri2c==1,]$s__Prevotella_copri)
#chagne
weightchg_blood.s__Prevotella_copri.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Prevotella_copri.uni)

weightchg_blood.s__Prevotella_copri.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Prevotella_copri.uni)
#chagne %
weightchgp_blood.s__Prevotella_copri.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                                + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Prevotella_copri.uni)

weightchgp_blood.s__Prevotella_copri.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                                +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Prevotella_copri.uni)

#6. s__Alistipes_putredinis
summary(meta_species_12stool$s__Alistipes_putredinis)
table(meta_species_12stool$s__Alistipes_putredinis)
graphics.off()
histogram(meta_species_12stool$s__Alistipes_putredinis)
meta_species_12stool$s__Alistipes_putredinis2c <- with(meta_species_12stool, ifelse(s__Alistipes_putredinis>3.1016 , 1, 0))
table(meta_species_12stool$s__Alistipes_putredinis2c)
summary(meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1,]$s__Alistipes_putredinis)
#chagne
weightchg_blood.s__Alistipes_putredinis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_putredinis.uni)
# all act 
weightchg_blood.s__Alistipes_putredinis.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_putredinis.uni)

# all act in each stratum
lmer.weightchg_blood.s__Alistipes_putredinis.below <- lmer ( weightchg_blood ~ paee_pamwk 
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.below)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.below)

lmer.weightchg_blood.s__Alistipes_putredinis.above <- lmer ( weightchg_blood ~ paee_pamwk
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.above)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.above)

# vig act 
weightchg_blood.s__Alistipes_putredinis.uni <- lmer( weightchg_blood ~ mets_vigpamwk + s__Alistipes_putredinis2c + mets_vigpamwk*s__Alistipes_putredinis2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_putredinis.uni)

# vig act in each stratum
lmer.weightchg_blood.s__Alistipes_putredinis.below <- lmer ( weightchg_blood ~ mets_vigpamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.below)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.below)

lmer.weightchg_blood.s__Alistipes_putredinis.above <- lmer ( weightchg_blood ~ mets_vigpamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.above)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.above)

# mod act 
weightchg_blood.s__Alistipes_putredinis.uni <- lmer( weightchg_blood ~ mets_modpamwk + s__Alistipes_putredinis2c + mets_modpamwk*s__Alistipes_putredinis2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_putredinis.uni)

# mod act in each stratum
lmer.weightchg_blood.s__Alistipes_putredinis.below <- lmer ( weightchg_blood ~ mets_modpamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.below)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.below)

lmer.weightchg_blood.s__Alistipes_putredinis.above <- lmer ( weightchg_blood ~ mets_modpamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.above)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.above)

# light act 
weightchg_blood.s__Alistipes_putredinis.uni <- lmer( weightchg_blood ~ mets_ltpamwk + s__Alistipes_putredinis2c + mets_ltpamwk*s__Alistipes_putredinis2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Alistipes_putredinis.uni)

# light act in each stratum
lmer.weightchg_blood.s__Alistipes_putredinis.below <- lmer ( weightchg_blood ~ mets_ltpamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==0, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.below)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.below)

lmer.weightchg_blood.s__Alistipes_putredinis.above <- lmer ( weightchg_blood ~ mets_ltpamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1, ])
summary(lmer.weightchg_blood.s__Alistipes_putredinis.above)
confint(lmer.weightchg_blood.s__Alistipes_putredinis.above)

#chagne %
weightchgp_blood.s__Alistipes_putredinis.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Alistipes_putredinis.uni)

weightchgp_blood.s__Alistipes_putredinis.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Alistipes_putredinis.uni)

# in each stratum
lmer.weightchgp_blood.s__Alistipes_putredinis.below <- lmer ( weightchgp_blood ~ paee_pamwk 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==0, ])
summary(lmer.weightchgp_blood.s__Alistipes_putredinis.below)
confint(lmer.weightchgp_blood.s__Alistipes_putredinis.below)

lmer.weightchgp_blood.s__Alistipes_putredinis.above <- lmer ( weightchgp_blood ~ paee_pamwk
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species_12stool[meta_species_12stool$s__Alistipes_putredinis2c==1, ])
summary(lmer.weightchgp_blood.s__Alistipes_putredinis.above)
confint(lmer.weightchgp_blood.s__Alistipes_putredinis.above)

#7. s__Eubacterium_siraeum
summary(meta_species_12stool$s__Eubacterium_siraeum)
table(meta_species_12stool$s__Eubacterium_siraeum)
graphics.off()
histogram(meta_species_12stool$s__Eubacterium_siraeum)
meta_species_12stool$s__Eubacterium_siraeum2c <- with(meta_species_12stool, ifelse(s__Eubacterium_siraeum> 0.2039, 1, 0))
table(meta_species_12stool$s__Eubacterium_siraeum2c)
summary(meta_species_12stool[meta_species_12stool$s__Eubacterium_siraeum2c==1,]$s__Eubacterium_siraeum)
#chagne
weightchg_blood.s__Eubacterium_siraeum.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Eubacterium_siraeum.uni)

weightchg_blood.s__Eubacterium_siraeum.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                                  +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Eubacterium_siraeum.uni)
#chagne %
weightchgp_blood.s__Eubacterium_siraeum.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                                   + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Eubacterium_siraeum.uni)

weightchgp_blood.s__Eubacterium_siraeum.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Eubacterium_siraeum.uni)

#8. s__Ruminococcus_bromii
summary(meta_species_12stool$s__Ruminococcus_bromii)
table(meta_species_12stool$s__Ruminococcus_bromii)
graphics.off()
histogram(meta_species_12stool$s__Ruminococcus_bromii)
meta_species_12stool$s__Ruminococcus_bromii2c <- with(meta_species_12stool, ifelse(s__Ruminococcus_bromii> 0.7604, 1, 0))
table(meta_species_12stool$s__Ruminococcus_bromii2c)
summary(meta_species_12stool[meta_species_12stool$s__Ruminococcus_bromii2c==1,]$s__Ruminococcus_bromii)
#chagne
weightchg_blood.s__Ruminococcus_bromii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Ruminococcus_bromii.uni)

weightchg_blood.s__Ruminococcus_bromii.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                                  +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Ruminococcus_bromii.uni)
#chagne %
weightchgp_blood.s__Ruminococcus_bromii.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                                   + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Ruminococcus_bromii.uni)

weightchgp_blood.s__Ruminococcus_bromii.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Ruminococcus_bromii.uni)

#9. s__Bacteroides_vulgatus
summary(meta_species_12stool$s__Bacteroides_vulgatus)
table(meta_species_12stool$s__Bacteroides_vulgatus)
graphics.off()
histogram(meta_species_12stool$s__Bacteroides_vulgatus)
meta_species_12stool$s__Bacteroides_vulgatus2c <- with(meta_species_12stool, ifelse(s__Bacteroides_vulgatus> 1.544, 1, 0))
table(meta_species_12stool$s__Bacteroides_vulgatus2c)
summary(meta_species_12stool[meta_species_12stool$s__Bacteroides_vulgatus2c==1,]$s__Bacteroides_vulgatus)
#chagne
weightchg_blood.s__Bacteroides_vulgatus.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Bacteroides_vulgatus.uni)

weightchg_blood.s__Bacteroides_vulgatus.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Bacteroides_vulgatus.uni)
#chagne %

weightchgp_blood.s__Bacteroides_vulgatus.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                                    + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Bacteroides_vulgatus.uni)

weightchgp_blood.s__Bacteroides_vulgatus.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Bacteroides_vulgatus.uni)

#10. s__Bacteroides_stercoris
summary(meta_species_12stool$s__Bacteroides_stercoris)
table(meta_species_12stool$s__Bacteroides_stercoris)
graphics.off()
histogram(meta_species_12stool$s__Bacteroides_stercoris)
meta_species_12stool$s__Bacteroides_stercoris2c <- with(meta_species_12stool, ifelse(s__Bacteroides_stercoris> 0.00672, 1, 0))
table(meta_species_12stool$s__Bacteroides_stercoris2c)
summary(meta_species_12stool[meta_species_12stool$s__Bacteroides_stercoris2c==1,]$s__Bacteroides_stercoris)
#chagne
weightchg_blood.s__Bacteroides_stercoris.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                                           + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Bacteroides_stercoris.uni)

weightchg_blood.s__Bacteroides_stercoris.uni <- lmer( weightchg_blood ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchg_blood.s__Bacteroides_stercoris.uni)
#chagne %
weightchgp_blood.s__Bacteroides_stercoris.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                                     + (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Bacteroides_stercoris.uni)

weightchgp_blood.s__Bacteroides_stercoris.uni <- lmer( weightchgp_blood ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species_12stool)
summary(weightchgp_blood.s__Bacteroides_stercoris.uni)
