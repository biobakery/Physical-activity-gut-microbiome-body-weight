
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity and the top 10 abundant species in relation to body weight change
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
meta_species$mlpa_paqlong <- meta_species$mpa_paqlong+meta_species$lpa_paqlong

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

####### bmi_dlw ########
#1. s__Eubacterium_rectale
summary(meta_species$s__Eubacterium_rectale)
table(meta_species$s__Eubacterium_rectale)
graphics.off()
histogram(meta_species$s__Eubacterium_rectale)
meta_species$s__Eubacterium_rectale2c <- with(meta_species, ifelse(s__Eubacterium_rectale> 6.616, 1, 0))
table(meta_species$s__Eubacterium_rectale2c)
summary(meta_species[meta_species$s__Eubacterium_rectale2c==1,]$s__Eubacterium_rectale)

bmi_dlw.s__Eubacterium_rectale.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                          + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Eubacterium_rectale.uni)

bmi_dlw.s__Eubacterium_rectale.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                                  +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Eubacterium_rectale.uni)

#2. s__Faecalibacterium_prausnitzii
summary(meta_species$s__Faecalibacterium_prausnitzii)
table(meta_species$s__Faecalibacterium_prausnitzii)
graphics.off()
histogram(meta_species$s__Faecalibacterium_prausnitzii)
meta_species$s__Faecalibacterium_prausnitzii2c <- with(meta_species, ifelse(s__Faecalibacterium_prausnitzii> 6.452, 1, 0))
table(meta_species$s__Faecalibacterium_prausnitzii2c)
summary(meta_species[meta_species$s__Faecalibacterium_prausnitzii2c==1,]$s__Faecalibacterium_prausnitzii)

bmi_dlw.s__Faecalibacterium_prausnitzii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                  + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Faecalibacterium_prausnitzii.uni)

bmi_dlw.s__Faecalibacterium_prausnitzii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Faecalibacterium_prausnitzii.uni)

#3. s__Subdoligranulum_unclassified
summary(meta_species$s__Subdoligranulum_unclassified)
table(meta_species$s__Subdoligranulum_unclassified)
graphics.off()
histogram(meta_species$s__Subdoligranulum_unclassified)
meta_species$s__Subdoligranulum_unclassified2c <- with(meta_species, ifelse(s__Subdoligranulum_unclassified> 5.234, 1, 0))
table(meta_species$s__Subdoligranulum_unclassified2c)
summary(meta_species[meta_species$s__Subdoligranulum_unclassified2c==1,]$s__Subdoligranulum_unclassified)

bmi_dlw.s__Subdoligranulum_unclassified.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Subdoligranulum_unclassified.uni)

bmi_dlw.s__Subdoligranulum_unclassified.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                           +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Subdoligranulum_unclassified.uni)

#4. s__Bacteroides_uniformis
summary(meta_species$s__Bacteroides_uniformis)
table(meta_species$s__Bacteroides_uniformis)
graphics.off()
histogram(meta_species$s__Bacteroides_uniformis)
meta_species$s__Bacteroides_uniformis2c <- with(meta_species, ifelse(s__Bacteroides_uniformis> 4.059, 1, 0))
table(meta_species$s__Bacteroides_uniformis2c)
summary(meta_species[meta_species$s__Bacteroides_uniformis2c==1,]$s__Bacteroides_uniformis)

bmi_dlw.s__Bacteroides_uniformis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Bacteroides_uniformis.uni)

bmi_dlw.s__Bacteroides_uniformis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Bacteroides_uniformis.uni)

#5. s__Prevotella_copri
summary(meta_species$s__Prevotella_copri)
table(meta_species$s__Prevotella_copri)
graphics.off()
histogram(meta_species$s__Prevotella_copri)
meta_species$s__Prevotella_copri2c <- with(meta_species, ifelse(s__Prevotella_copri> 0, 1, 0))
table(meta_species$s__Prevotella_copri2c)
summary(meta_species[meta_species$s__Prevotella_copri2c==1,]$s__Prevotella_copri)

bmi_dlw.s__Prevotella_copri.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Prevotella_copri.uni)

bmi_dlw.s__Prevotella_copri.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                               +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Prevotella_copri.uni)

#6. s__Alistipes_putredinis
summary(meta_species$s__Alistipes_putredinis)
table(meta_species$s__Alistipes_putredinis)
graphics.off()
histogram(meta_species$s__Alistipes_putredinis)
meta_species$s__Alistipes_putredinis2c <- with(meta_species, ifelse(s__Alistipes_putredinis> 3.0033, 1, 0))
table(meta_species$s__Alistipes_putredinis2c)
summary(meta_species[meta_species$s__Alistipes_putredinis2c==1,]$s__Alistipes_putredinis)

bmi_dlw.s__Alistipes_putredinis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_putredinis.uni)

bmi_dlw.s__Alistipes_putredinis.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_putredinis.uni)

#all act in each stratum
lmer.bmi_dlw.s__Alistipes_putredinis.below <- lmer ( bmi_dlw ~ paee_pamwk 
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.below)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.below)

lmer.bmi_dlw.s__Alistipes_putredinis.above <- lmer ( bmi_dlw ~ paee_pamwk
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.above)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.above)

#vig act
bmi_dlw.s__Alistipes_putredinis.uni <- lmer( bmi_dlw ~ mets_vigpamwk + s__Alistipes_putredinis2c + mets_vigpamwk*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_putredinis.uni)

#vig act in each stratum
lmer.bmi_dlw.s__Alistipes_putredinis.below <- lmer ( bmi_dlw ~ mets_vigpamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.below)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.below)

lmer.bmi_dlw.s__Alistipes_putredinis.above <- lmer ( bmi_dlw ~ mets_vigpamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.above)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.above)

#mod act 
bmi_dlw.s__Alistipes_putredinis.uni <- lmer( bmi_dlw ~ mets_modpamwk + s__Alistipes_putredinis2c + mets_modpamwk*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_putredinis.uni)

#mod act in each stratum
lmer.bmi_dlw.s__Alistipes_putredinis.below <- lmer ( bmi_dlw ~ mets_modpamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.below)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.below)

lmer.bmi_dlw.s__Alistipes_putredinis.above <- lmer ( bmi_dlw ~ mets_modpamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.above)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.above)

#light act
bmi_dlw.s__Alistipes_putredinis.uni <- lmer( bmi_dlw ~ mets_ltpamwk + s__Alistipes_putredinis2c + mets_ltpamwk*s__Alistipes_putredinis2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Alistipes_putredinis.uni)

#light act in each stratum
lmer.bmi_dlw.s__Alistipes_putredinis.below <- lmer ( bmi_dlw ~ mets_ltpamwk 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.below)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.below)

lmer.bmi_dlw.s__Alistipes_putredinis.above <- lmer ( bmi_dlw ~ mets_ltpamwk
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.bmi_dlw.s__Alistipes_putredinis.above)
confint(lmer.bmi_dlw.s__Alistipes_putredinis.above)

#7. s__Eubacterium_siraeum
summary(meta_species$s__Eubacterium_siraeum)
table(meta_species$s__Eubacterium_siraeum)
graphics.off()
histogram(meta_species$s__Eubacterium_siraeum)
meta_species$s__Eubacterium_siraeum2c <- with(meta_species, ifelse(s__Eubacterium_siraeum> 0.2456, 1, 0))
table(meta_species$s__Eubacterium_siraeum2c)
summary(meta_species[meta_species$s__Eubacterium_siraeum2c==1,]$s__Eubacterium_siraeum)

bmi_dlw.s__Eubacterium_siraeum.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Eubacterium_siraeum.uni)

bmi_dlw.s__Eubacterium_siraeum.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                                  +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Eubacterium_siraeum.uni)

#8. s__Ruminococcus_bromii
summary(meta_species$s__Ruminococcus_bromii)
table(meta_species$s__Ruminococcus_bromii)
graphics.off()
histogram(meta_species$s__Ruminococcus_bromii)
meta_species$s__Ruminococcus_bromii2c <- with(meta_species, ifelse(s__Ruminococcus_bromii> 0.6912, 1, 0))
table(meta_species$s__Ruminococcus_bromii2c)
summary(meta_species[meta_species$s__Ruminococcus_bromii2c==1,]$s__Ruminococcus_bromii)

bmi_dlw.s__Ruminococcus_bromii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Ruminococcus_bromii.uni)

bmi_dlw.s__Ruminococcus_bromii.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                                  +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Ruminococcus_bromii.uni)

#9. s__Bacteroides_vulgatus
summary(meta_species$s__Bacteroides_vulgatus)
table(meta_species$s__Bacteroides_vulgatus)
graphics.off()
histogram(meta_species$s__Bacteroides_vulgatus)
meta_species$s__Bacteroides_vulgatus2c <- with(meta_species, ifelse(s__Bacteroides_vulgatus> 1.5450, 1, 0))
table(meta_species$s__Bacteroides_vulgatus2c)
summary(meta_species[meta_species$s__Bacteroides_vulgatus2c==1,]$s__Bacteroides_vulgatus)

bmi_dlw.s__Bacteroides_vulgatus.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Bacteroides_vulgatus.uni)

bmi_dlw.s__Bacteroides_vulgatus.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                                   +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Bacteroides_vulgatus.uni)

#10. s__Bacteroides_stercoris
summary(meta_species$s__Bacteroides_stercoris)
table(meta_species$s__Bacteroides_stercoris)
graphics.off()
histogram(meta_species$s__Bacteroides_stercoris)
meta_species$s__Bacteroides_stercoris2c <- with(meta_species, ifelse(s__Bacteroides_stercoris> 0.02914, 1, 0))
table(meta_species$s__Bacteroides_stercoris2c)
summary(meta_species[meta_species$s__Bacteroides_stercoris2c==1,]$s__Bacteroides_stercoris)

bmi_dlw.s__Bacteroides_stercoris.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Bacteroides_stercoris.uni)

bmi_dlw.s__Bacteroides_stercoris.uni <- lmer( bmi_dlw ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.s__Bacteroides_stercoris.uni)

####### pfat_dlw ########

#1. s__Eubacterium_rectale
pfat_dlw.s__Eubacterium_rectale.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                                  + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Eubacterium_rectale.uni)
pfat_dlw.s__Eubacterium_rectale.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Eubacterium_rectale.uni)

#2. s__Faecalibacterium_prausnitzii
pfat_dlw.s__Faecalibacterium_prausnitzii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Faecalibacterium_prausnitzii.uni)

pfat_dlw.s__Faecalibacterium_prausnitzii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Faecalibacterium_prausnitzii.uni)

#3. s__Subdoligranulum_unclassified
pfat_dlw.s__Subdoligranulum_unclassified.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                           + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Subdoligranulum_unclassified.uni)

pfat_dlw.s__Subdoligranulum_unclassified.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Subdoligranulum_unclassified.uni)

#4. s__Bacteroides_uniformis
pfat_dlw.s__Bacteroides_uniformis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                                    + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Bacteroides_uniformis.uni)

pfat_dlw.s__Bacteroides_uniformis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Bacteroides_uniformis.uni)

#5. s__Prevotella_copri

pfat_dlw.s__Prevotella_copri.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                               + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Prevotella_copri.uni)

pfat_dlw.s__Prevotella_copri.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Prevotella_copri.uni)

#6. s__Alistipes_putredinis
pfat_dlw.s__Alistipes_putredinis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                                   + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_putredinis.uni)
#all act
pfat_dlw.s__Alistipes_putredinis.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_putredinis.uni)

# all act in each stratum
lmer.pfat_dlw.s__Alistipes_putredinis.below <- lmer ( pfat_dlw ~ paee_pamwk 
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.below)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.below)

lmer.pfat_dlw.s__Alistipes_putredinis.above <- lmer ( pfat_dlw ~ paee_pamwk 
                                                          +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.above)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.above)

#vig act
pfat_dlw.s__Alistipes_putredinis.uni <- lmer( pfat_dlw ~ mets_vigpamwk + s__Alistipes_putredinis2c + mets_vigpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_putredinis.uni)

# vig act in each stratum
lmer.pfat_dlw.s__Alistipes_putredinis.below <- lmer ( pfat_dlw ~ mets_vigpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.below)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.below)

lmer.pfat_dlw.s__Alistipes_putredinis.above <- lmer ( pfat_dlw ~ mets_vigpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.above)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.above)

#mod act
pfat_dlw.s__Alistipes_putredinis.uni <- lmer( pfat_dlw ~ mets_modpamwk + s__Alistipes_putredinis2c + mets_modpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_putredinis.uni)

# mod act in each stratum
lmer.pfat_dlw.s__Alistipes_putredinis.below <- lmer ( pfat_dlw ~ mets_modpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.below)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.below)

lmer.pfat_dlw.s__Alistipes_putredinis.above <- lmer ( pfat_dlw ~ mets_modpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.above)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.above)

#light act
pfat_dlw.s__Alistipes_putredinis.uni <- lmer( pfat_dlw ~ mets_ltpamwk + s__Alistipes_putredinis2c + mets_ltpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Alistipes_putredinis.uni)

# light act in each stratum
lmer.pfat_dlw.s__Alistipes_putredinis.below <- lmer ( pfat_dlw ~ mets_ltpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.below)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.below)

lmer.pfat_dlw.s__Alistipes_putredinis.above <- lmer ( pfat_dlw ~ mets_ltpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.pfat_dlw.s__Alistipes_putredinis.above)
confint(lmer.pfat_dlw.s__Alistipes_putredinis.above)

#7. s__Eubacterium_siraeum
pfat_dlw.s__Eubacterium_siraeum.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                                  + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Eubacterium_siraeum.uni)

pfat_dlw.s__Eubacterium_siraeum.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Eubacterium_siraeum.uni)

#8. s__Ruminococcus_bromii
pfat_dlw.s__Ruminococcus_bromii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                                  + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Ruminococcus_bromii.uni)

pfat_dlw.s__Ruminococcus_bromii.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Ruminococcus_bromii.uni)

#9. s__Bacteroides_vulgatus

pfat_dlw.s__Bacteroides_vulgatus.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                                   + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Bacteroides_vulgatus.uni)

pfat_dlw.s__Bacteroides_vulgatus.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Bacteroides_vulgatus.uni)

#10. s__Bacteroides_stercoris
pfat_dlw.s__Bacteroides_stercoris.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                                    + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Bacteroides_stercoris.uni)

pfat_dlw.s__Bacteroides_stercoris.uni <- lmer( pfat_dlw ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.s__Bacteroides_stercoris.uni)

####### wtchgsto21 ########

#1. s__Eubacterium_rectale
wtchg6521.s__Eubacterium_rectale.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Eubacterium_rectale2c + act_paqlong*s__Eubacterium_rectale2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Eubacterium_rectale.uni)

wtchg6521.s__Eubacterium_rectale.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Eubacterium_rectale2c + act_paqlong*s__Eubacterium_rectale2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Eubacterium_rectale.uni)

#2. s__Faecalibacterium_prausnitzii
wtchg6521.s__Faecalibacterium_prausnitzii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Faecalibacterium_prausnitzii2c + act_paqlong*s__Faecalibacterium_prausnitzii2c 
                                                     + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Faecalibacterium_prausnitzii.uni)

wtchg6521.s__Faecalibacterium_prausnitzii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Faecalibacterium_prausnitzii2c + act_paqlong*s__Faecalibacterium_prausnitzii2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Faecalibacterium_prausnitzii.uni)

#3. s__Subdoligranulum_unclassified
wtchg6521.s__Subdoligranulum_unclassified.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Subdoligranulum_unclassified2c + act_paqlong*s__Subdoligranulum_unclassified2c 
                                                     + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Subdoligranulum_unclassified.uni)

wtchg6521.s__Subdoligranulum_unclassified.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Subdoligranulum_unclassified2c + act_paqlong*s__Subdoligranulum_unclassified2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Subdoligranulum_unclassified.uni)

#4. s__Bacteroides_uniformis
wtchg6521.s__Bacteroides_uniformis.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Bacteroides_uniformis2c + act_paqlong*s__Bacteroides_uniformis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Bacteroides_uniformis.uni)

wtchg6521.s__Bacteroides_uniformis.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Bacteroides_uniformis2c + act_paqlong*s__Bacteroides_uniformis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Bacteroides_uniformis.uni)

#5. s__Prevotella_copri
wtchg6521.s__Prevotella_copri.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Prevotella_copri2c + act_paqlong*s__Prevotella_copri2c 
                                         + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Prevotella_copri.uni)

wtchg6521.s__Prevotella_copri.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Prevotella_copri2c + act_paqlong*s__Prevotella_copri2c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Prevotella_copri.uni)

#6. s__Alistipes_putredinis
wtchg6521.s__Alistipes_putredinis.uni <- lmer( wtchgsto21 ~ act_paqlong_sec + s__Alistipes_putredinis2c + act_paqlong_sec*s__Alistipes_putredinis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Alistipes_putredinis.uni)

wtchg6521.s__Alistipes_putredinis.uni <- lmer( wtchgsto21 ~ act_paqlong_sec + s__Alistipes_putredinis2c + act_paqlong_sec*s__Alistipes_putredinis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Alistipes_putredinis.uni)

# all pa in each stratum
lmer.wtchg6521.s__Alistipes_putredinis.below <- lmer ( wtchgsto21 ~ act_paqlong_sec 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.below)
confint(lmer.wtchg6521.s__Alistipes_putredinis.below)

lmer.wtchg6521.s__Alistipes_putredinis.above <- lmer ( wtchgsto21 ~ act_paqlong_sec 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.above)
confint(lmer.wtchg6521.s__Alistipes_putredinis.above)

# vigorous pa 
wtchg6521.s__Alistipes_putredinis.uni <- lmer( wtchgsto21 ~ vpa_paqlong + s__Alistipes_putredinis2c + vpa_paqlong*s__Alistipes_putredinis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Alistipes_putredinis.uni)

# vigorous pa in each stratum
lmer.wtchg6521.s__Alistipes_putredinis.below <- lmer ( wtchgsto21 ~ vpa_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.below)
confint(lmer.wtchg6521.s__Alistipes_putredinis.below)

lmer.wtchg6521.s__Alistipes_putredinis.above <- lmer ( wtchgsto21 ~ vpa_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.above)
confint(lmer.wtchg6521.s__Alistipes_putredinis.above)

# moderate pa 
wtchg6521.s__Alistipes_putredinis.uni <- lmer( wtchgsto21 ~ mlpa_paqlong + s__Alistipes_putredinis2c + mlpa_paqlong*s__Alistipes_putredinis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Alistipes_putredinis.uni)

# moderate pa in each stratum
lmer.wtchg6521.s__Alistipes_putredinis.below <- lmer ( wtchgsto21 ~ mlpa_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.below)
confint(lmer.wtchg6521.s__Alistipes_putredinis.below)

lmer.wtchg6521.s__Alistipes_putredinis.above <- lmer ( wtchgsto21 ~ mlpa_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.above)
confint(lmer.wtchg6521.s__Alistipes_putredinis.above)

# light pa 
wtchg6521.s__Alistipes_putredinis.uni <- lmer( wtchgsto21 ~ lpa_paqlong + s__Alistipes_putredinis2c + lpa_paqlong*s__Alistipes_putredinis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Alistipes_putredinis.uni)

# light pa in each stratum
lmer.wtchg6521.s__Alistipes_putredinis.below <- lmer ( wtchgsto21 ~ lpa_paqlong+vpa_paqlong+mpa_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.below)
confint(lmer.wtchg6521.s__Alistipes_putredinis.below)

lmer.wtchg6521.s__Alistipes_putredinis.above <- lmer ( wtchgsto21 ~ lpa_paqlong+vpa_paqlong+mpa_paqlong 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.wtchg6521.s__Alistipes_putredinis.above)
confint(lmer.wtchg6521.s__Alistipes_putredinis.above)

#7. s__Eubacterium_siraeum
wtchg6521.s__Eubacterium_siraeum.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Eubacterium_siraeum2c + act_paqlong*s__Eubacterium_siraeum2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Eubacterium_siraeum.uni)

wtchg6521.s__Eubacterium_siraeum.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Eubacterium_siraeum2c + act_paqlong*s__Eubacterium_siraeum2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Eubacterium_siraeum.uni)

#8. s__Ruminococcus_bromii
wtchg6521.s__Ruminococcus_bromii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Ruminococcus_bromii2c + act_paqlong*s__Ruminococcus_bromii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Ruminococcus_bromii.uni)

wtchg6521.s__Ruminococcus_bromii.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Ruminococcus_bromii2c + act_paqlong*s__Ruminococcus_bromii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Ruminococcus_bromii.uni)

#9. s__Bacteroides_vulgatus
wtchg6521.s__Bacteroides_vulgatus.uni <- lmer( wtchgsto21 ~ act_paqlong+ s__Bacteroides_vulgatus2c + act_paqlong*s__Bacteroides_vulgatus2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Bacteroides_vulgatus.uni)

wtchg6521.s__Bacteroides_vulgatus.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Bacteroides_vulgatus2c + act_paqlong*s__Bacteroides_vulgatus2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Bacteroides_vulgatus.uni)

#10. s__Bacteroides_stercoris
wtchg6521.s__Bacteroides_stercoris.uni <- lmer( wtchgsto21 ~ act_paqlong + s__Bacteroides_stercoris2c + act_paqlong*s__Bacteroides_stercoris2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Bacteroides_stercoris.uni)

wtchg6521.s__Bacteroides_stercoris.uni <- lmer( wtchgsto21~ act_paqlong + s__Bacteroides_stercoris2c + act_paqlong*s__Bacteroides_stercoris2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(wtchg6521.s__Bacteroides_stercoris.uni)

####### lg2hba1c ########

#1. s__Eubacterium_rectale
lg2hba1c.s__Eubacterium_rectale.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Eubacterium_rectale.uni)

lg2hba1c.s__Eubacterium_rectale.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Eubacterium_rectale.uni)

#2. s__Faecalibacterium_prausnitzii
lg2hba1c.s__Faecalibacterium_prausnitzii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                     + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Faecalibacterium_prausnitzii.uni)

lg2hba1c.s__Faecalibacterium_prausnitzii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Faecalibacterium_prausnitzii.uni)

#3. s__Subdoligranulum_unclassified
lg2hba1c.s__Subdoligranulum_unclassified.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                     + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Subdoligranulum_unclassified.uni)

lg2hba1c.s__Subdoligranulum_unclassified.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Subdoligranulum_unclassified.uni)

#4. s__Bacteroides_uniformis
lg2hba1c.s__Bacteroides_uniformis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Bacteroides_uniformis.uni)

lg2hba1c.s__Bacteroides_uniformis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Bacteroides_uniformis.uni)

#5. s__Prevotella_copri
lg2hba1c.s__Prevotella_copri.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                         + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Prevotella_copri.uni)

lg2hba1c.s__Prevotella_copri.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Prevotella_copri.uni)

#6. s__Alistipes_putredinis
lg2hba1c.s__Alistipes_putredinis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_putredinis.uni)

lg2hba1c.s__Alistipes_putredinis.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_putredinis.uni)

# all act in each stratum
lmer.lg2hba1c.s__Alistipes_putredinis.below <- lmer ( lg2hba1c ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.below)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.below)

lmer.lg2hba1c.s__Alistipes_putredinis.above <- lmer ( lg2hba1c ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.above)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.above)

# vig act
lg2hba1c.s__Alistipes_putredinis.uni <- lmer( lg2hba1c ~ mets_vigpamwk + s__Alistipes_putredinis2c + mets_vigpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_putredinis.uni)

# vig act in each stratum
lmer.lg2hba1c.s__Alistipes_putredinis.below <- lmer ( lg2hba1c ~ mets_vigpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.below)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.below)

lmer.lg2hba1c.s__Alistipes_putredinis.above <- lmer ( lg2hba1c ~ mets_vigpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.above)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.above)

# mod act 
lg2hba1c.s__Alistipes_putredinis.uni <- lmer( lg2hba1c ~ mets_modpamwk + s__Alistipes_putredinis2c + mets_modpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_putredinis.uni)

# mod act in each stratum
lmer.lg2hba1c.s__Alistipes_putredinis.below <- lmer ( lg2hba1c ~ mets_modpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.below)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.below)

lmer.lg2hba1c.s__Alistipes_putredinis.above <- lmer ( lg2hba1c ~ mets_modpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.above)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.above)

# light act 
lg2hba1c.s__Alistipes_putredinis.uni <- lmer( lg2hba1c ~ mets_ltpamwk + s__Alistipes_putredinis2c + mets_ltpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Alistipes_putredinis.uni)

# light act in each stratum
lmer.lg2hba1c.s__Alistipes_putredinis.below <- lmer ( lg2hba1c ~ mets_ltpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.below)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.below)

lmer.lg2hba1c.s__Alistipes_putredinis.above <- lmer ( lg2hba1c ~ mets_ltpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2hba1c.s__Alistipes_putredinis.above)
confint(lmer.lg2hba1c.s__Alistipes_putredinis.above)

#7. s__Eubacterium_siraeum
lg2hba1c.s__Eubacterium_siraeum.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Eubacterium_siraeum.uni)

lg2hba1c.s__Eubacterium_siraeum.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Eubacterium_siraeum.uni)

#8. s__Ruminococcus_bromii
lg2hba1c.s__Ruminococcus_bromii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Ruminococcus_bromii.uni)

lg2hba1c.s__Ruminococcus_bromii.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Ruminococcus_bromii.uni)

#9. s__Bacteroides_vulgatus
lg2hba1c.s__Bacteroides_vulgatus.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Bacteroides_vulgatus.uni)

lg2hba1c.s__Bacteroides_vulgatus.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Bacteroides_vulgatus.uni)

#10. s__Bacteroides_stercoris
lg2hba1c.s__Bacteroides_stercoris.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Bacteroides_stercoris.uni)

lg2hba1c.s__Bacteroides_stercoris.uni <- lmer( lg2hba1c ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2hba1c.s__Bacteroides_stercoris.uni)

####### lg2crp ########

#1. s__Eubacterium_rectale
lg2crp.s__Eubacterium_rectale.uni <- lmer( lg2crp ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Eubacterium_rectale.uni)

lg2crp.s__Eubacterium_rectale.uni <- lmer( lg2crp ~ paee_pamwk + s__Eubacterium_rectale2c + paee_pamwk*s__Eubacterium_rectale2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Eubacterium_rectale.uni)

#2. s__Faecalibacterium_prausnitzii
lg2crp.s__Faecalibacterium_prausnitzii.uni <- lmer( lg2crp ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                     + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Faecalibacterium_prausnitzii.uni)

lg2crp.s__Faecalibacterium_prausnitzii.uni <- lmer( lg2crp ~ paee_pamwk + s__Faecalibacterium_prausnitzii2c + paee_pamwk*s__Faecalibacterium_prausnitzii2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Faecalibacterium_prausnitzii.uni)

#3. s__Subdoligranulum_unclassified
lg2crp.s__Subdoligranulum_unclassified.uni <- lmer( lg2crp ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                     + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Subdoligranulum_unclassified.uni)

lg2crp.s__Subdoligranulum_unclassified.uni <- lmer( lg2crp ~ paee_pamwk + s__Subdoligranulum_unclassified2c + paee_pamwk*s__Subdoligranulum_unclassified2c 
                                                     +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Subdoligranulum_unclassified.uni)

#4. s__Bacteroides_uniformis
lg2crp.s__Bacteroides_uniformis.uni <- lmer( lg2crp ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Bacteroides_uniformis.uni)

lg2crp.s__Bacteroides_uniformis.uni <- lmer( lg2crp ~ paee_pamwk + s__Bacteroides_uniformis2c + paee_pamwk*s__Bacteroides_uniformis2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Bacteroides_uniformis.uni)

#5. s__Prevotella_copri
lg2crp.s__Prevotella_copri.uni <- lmer( lg2crp ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                         + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Prevotella_copri.uni)

lg2crp.s__Prevotella_copri.uni <- lmer( lg2crp ~ paee_pamwk + s__Prevotella_copri2c + paee_pamwk*s__Prevotella_copri2c 
                                         +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Prevotella_copri.uni)

#6. s__Alistipes_putredinis
# all act 
lg2crp.s__Alistipes_putredinis.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_putredinis.uni)

lg2crp.s__Alistipes_putredinis.uni <- lmer( lg2crp ~ paee_pamwk + s__Alistipes_putredinis2c + paee_pamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_putredinis.uni)

# all act in each stratum
lmer.lg2crp.s__Alistipes_putredinis.below <- lmer ( lg2crp ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.below)
confint(lmer.lg2crp.s__Alistipes_putredinis.below)

lmer.lg2crp.s__Alistipes_putredinis.above <- lmer ( lg2crp ~ paee_pamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.above)
confint(lmer.lg2crp.s__Alistipes_putredinis.above)

# vig act 
lg2crp.s__Alistipes_putredinis.uni <- lmer( lg2crp ~ mets_vigpamwk + s__Alistipes_putredinis2c + mets_vigpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_putredinis.uni)

# vig act in each stratum
lmer.lg2crp.s__Alistipes_putredinis.below <- lmer ( lg2crp ~ mets_vigpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.below)
confint(lmer.lg2crp.s__Alistipes_putredinis.below)

lmer.lg2crp.s__Alistipes_putredinis.above <- lmer ( lg2crp ~ mets_vigpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.above)
confint(lmer.lg2crp.s__Alistipes_putredinis.above)

# mod act 
lg2crp.s__Alistipes_putredinis.uni <- lmer( lg2crp ~ mets_modpamwk + s__Alistipes_putredinis2c + mets_modpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_putredinis.uni)

# mod act in each stratum
lmer.lg2crp.s__Alistipes_putredinis.below <- lmer ( lg2crp ~ mets_modpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.below)
confint(lmer.lg2crp.s__Alistipes_putredinis.below)

lmer.lg2crp.s__Alistipes_putredinis.above <- lmer ( lg2crp ~ mets_modpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.above)
confint(lmer.lg2crp.s__Alistipes_putredinis.above)

# light act 
lg2crp.s__Alistipes_putredinis.uni <- lmer( lg2crp ~ mets_ltpamwk + s__Alistipes_putredinis2c + mets_ltpamwk*s__Alistipes_putredinis2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Alistipes_putredinis.uni)

# light act in each stratum
lmer.lg2crp.s__Alistipes_putredinis.below <- lmer ( lg2crp ~ mets_ltpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.below)
confint(lmer.lg2crp.s__Alistipes_putredinis.below)

lmer.lg2crp.s__Alistipes_putredinis.above <- lmer ( lg2crp ~ mets_ltpamwk 
                                                    +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(lmer.lg2crp.s__Alistipes_putredinis.above)
confint(lmer.lg2crp.s__Alistipes_putredinis.above)

#7. s__Eubacterium_siraeum
lg2crp.s__Eubacterium_siraeum.uni <- lmer( lg2crp ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Eubacterium_siraeum.uni)

lg2crp.s__Eubacterium_siraeum.uni <- lmer( lg2crp ~ paee_pamwk + s__Eubacterium_siraeum2c + paee_pamwk*s__Eubacterium_siraeum2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Eubacterium_siraeum.uni)

#8. s__Ruminococcus_bromii
lg2crp.s__Ruminococcus_bromii.uni <- lmer( lg2crp ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                            + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Ruminococcus_bromii.uni)

lg2crp.s__Ruminococcus_bromii.uni <- lmer( lg2crp ~ paee_pamwk + s__Ruminococcus_bromii2c + paee_pamwk*s__Ruminococcus_bromii2c 
                                            +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Ruminococcus_bromii.uni)

#9. s__Bacteroides_vulgatus
lg2crp.s__Bacteroides_vulgatus.uni <- lmer( lg2crp ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                             + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Bacteroides_vulgatus.uni)

lg2crp.s__Bacteroides_vulgatus.uni <- lmer( lg2crp ~ paee_pamwk + s__Bacteroides_vulgatus2c + paee_pamwk*s__Bacteroides_vulgatus2c 
                                             +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Bacteroides_vulgatus.uni)

#10. s__Bacteroides_stercoris
lg2crp.s__Bacteroides_stercoris.uni <- lmer( lg2crp ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                              + (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Bacteroides_stercoris.uni)

lg2crp.s__Bacteroides_stercoris.uni <- lmer( lg2crp ~ paee_pamwk + s__Bacteroides_stercoris2c + paee_pamwk*s__Bacteroides_stercoris2c 
                                              +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_species)
summary(lg2crp.s__Bacteroides_stercoris.uni)
