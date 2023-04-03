
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity and PCo1 and PCo2 in relation to body weight change
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
library(lmerTest)
library(lme4)

meta_species<-read.csv(file="./data_generated/meta_species.csv", header = TRUE)
meta_species$paee_pamwk <- meta_species$paee_pam*7

#bmi_dlw
bmi_dlw.paee_pamX1.uni <- lmer( bmi_dlw ~ paee_pam + X1 + paee_pam*X1 + (1 | SubjectID) , data=meta_species)
summary(bmi_dlw.paee_pamX1.uni)
bmi_dlw.paee_pamX1.uni <- lmer( bmi_dlw ~ paee_pam + X1 + paee_pam*X1 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species)
summary(bmi_dlw.paee_pamX1.uni)

#pfat_dlw
pfat_dlw.paee_pamX1.uni <- lmer( pfat_dlw ~ paee_pam + X1 + paee_pam*X1 + (1 | SubjectID) , data=meta_species)
summary(pfat_dlw.paee_pamX1.uni)
pfat_dlw.paee_pamX1.uni <- lmer( pfat_dlw ~ paee_pam + X1 + paee_pam*X1 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species)
summary(pfat_dlw.paee_pamX1.uni)

#wtchg_paq
wtchgp_paq.paee_pamwkX1.uni <- lmer( wtchg_paq ~ paee_pamwk + X2 + paee_pamwk*X2 + (1 | SubjectID) , data=meta_species_12stool)
summary(wtchgp_paq.paee_pamwkX1.uni)
wtchgp_paq.paee_pamwkX1.uni <- lmer( wtchg_paq ~ paee_pamwk + X2 + paee_pamwk*X2 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species_12stool)
summary(wtchgp_paq.paee_pamwkX1.uni)

#wtchgsto21
wtchgsto21.act_paqlongX1.uni <- lmer( wtchgsto21 ~ act_paqlong + X2 + act_paqlong*X2 + (1 | SubjectID) , data=meta_species)
summary(wtchgsto21.act_paqlongX1.uni)
wtchgsto21.act_paqlongX1.uni <- lmer( wtchgsto21 ~ act_paqlong + X2 + act_paqlong*X2 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_species)
summary(wtchgsto21.act_paqlongX1.uni)

#lg2hba1c
lg2hdl.paee_pamX1.uni <- lmer( lg2hba1c ~ paee_pamwk + X2 + paee_pamwk*X2 + (1 | SubjectID) , data=meta_species)
summary(lg2hdl.paee_pamX1.uni)
lg2hdl.paee_pamX1.uni <- lmer( lg2hba1c ~ paee_pamwk + X2 + paee_pamwk*X2 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species)
summary(lg2hdl.paee_pamX1.uni)

#lg2crp
lg2crp.paee_pamX1.uni <- lmer( lg2crp ~ paee_pamwk + X2 + paee_pamwk*X2 + (1 | SubjectID) , data=meta_species)
summary(lg2crp.paee_pamX1.uni)
lg2crp.paee_pamX1.uni <- lmer( lg2crp ~ paee_pamwk + X2 + paee_pamwk*X2 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species)
summary(lg2crp.paee_pamX1.uni)

#lg2tg
lg2tg.paee_pamX1.uni <- lmer( lg2tg ~ paee_pamwk + X2 + paee_pamwk*X2 + (1 | SubjectID) , data=meta_species)
summary(lg2tg.paee_pamX1.uni)
lg2tg.paee_pamX1.uni <- lmer( lg2tg ~ paee_pamwk + X2 + paee_pamwk*X2 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species)
summary(lg2tg.paee_pamX1.uni)

#lg2tc
lg2tc.paee_pamX1.uni <- lmer( lg2tc ~ paee_pamwk + X1 + paee_pamwk*X1 + (1 | SubjectID) , data=meta_species)
summary(lg2tc.paee_pamX1.uni)
lg2tc.paee_pamX1.uni <- lmer( lg2tc ~ paee_pamwk + X1 + paee_pamwk*X1 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species)
summary(lg2tc.paee_pamX1.uni)

#lg2hdl 
lg2hdl.paee_pamX1.uni <- lmer( lg2hdl ~ paee_pamwk + X2 + paee_pamwk*X2 + (1 | SubjectID) , data=meta_species)
summary(lg2hdl.paee_pamX1.uni)
lg2hdl.paee_pamX1.uni <- lmer( lg2hdl ~ paee_pamwk + X2 + paee_pamwk*X2 +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1|SubjectID) , data=meta_species)
summary(lg2hdl.paee_pamX1.uni)
