
##############################################################################################################################################
# 1) Purposes: conduct MaAsLin regression to calculate the associations of the measures of physical activity and adiposity with abundances of per microbial species
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

library(Maaslin2)
library(ggplot2)
library(grid)
library(ggpubr)
library(tidyverse)
library(gridExtra)

# read in metadata
meta_med<-read.csv(file="./data_generated/meta_species.csv")
meta_med$mvpa_paqlong <- meta_med$vpa_paqlong + meta_med$mpa_paqlong

# read in taxonomy data
tax_rpk_name <-   read.table(file = "./data_generated/bugs_dna_929_unFilt.tsv",
                             sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)

# only keep species-level features
tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))
rownames(tax_rpk_species)<-tax_rpk_species$species
tax_rpk_species<-tax_rpk_species[,-c(1:8)]

rownames(meta_med)<-meta_med$mid
meta_med$age_fecal<-as.numeric(as.character(meta_med$age_fecal))
meta_med$calor122cn<-as.numeric(as.character(meta_med$calor122cn))

# maaslin regressions
maaslin_diet<-function(dir, exposure){
  Maaslin2(tax_rpk_species, 
           meta_med, 
           dir,
           transform="AST",
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="SubjectID", 
           fixed_effects = c(exposure, "age_fecal","ahei","calor122cn", "probio_2mo_qu", "stool_type", "smk", "ant_12mo_qu"))
}

maaslin_diet(dir="./data_generated/maaslin_results/paee_pam_taxonomy",  exposure = "paee_pam")
maaslin_diet(dir="./data_generated/maaslin_results/mets_vigpam_taxonomy",  exposure = "mets_vigpam")
maaslin_diet(dir="./data_generated/maaslin_results/mets_modpam_taxonomy",  exposure = "mets_modpam")
maaslin_diet(dir="./data_generated/maaslin_results/mets_ltpam_taxonomy",  exposure = "mets_ltpam")
maaslin_diet(dir="./data_generated/maaslin_results/act_paqlong_taxonomy",  exposure = "act_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/vpa_paqlong_taxonomy",  exposure = "vpa_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/mpa_paqlong_taxonomy",  exposure = "mpa_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/lpa_paqlong_taxonomy",  exposure = "lpa_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/bmi_dlw_taxonomy",  exposure = "bmi_dlw")
maaslin_diet(dir="./data_generated/maaslin_results/pfat_dlw_taxonomy",  exposure = "pfat_dlw")
maaslin_diet(dir="./data_generated/maaslin_results/weightchg_blood_taxonomy",  exposure = "weightchg_blood")
maaslin_diet(dir="./data_generated/maaslin_results/wtchgsto21_taxonomy",  exposure = "wtchgsto21")
maaslin_diet(dir="./data_generated/maaslin_results/lgcrp_taxonomy",  exposure = "lgcrp")
maaslin_diet(dir="./data_generated/maaslin_results/lghba1c_taxonomy",  exposure = "lghba1c")


