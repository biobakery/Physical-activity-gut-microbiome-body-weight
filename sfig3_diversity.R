
##############################################################################################################################################
# 1) Purposes: calculated beta-diversity index (Shannon) and compared Shannon index across physical activity levels
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

library(ggplot2)
library(scales)
library(tidyverse)
library(vegan)
library(lmerTest)
library(lme4)

# read in taxonomy data
tax_rpk_name <-   read.table(    file = './bugs_dna_929_unFilt.tsv',
                                 sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)
# only keep species-level features
tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))
rownames(tax_rpk_species)<-tax_rpk_species$species
tax_rpk_species<-tax_rpk_species[,-c(1:8)]
ttax_rpk <- as.data.frame(t(tax_rpk_species))
ttax_rpk$shannon <- diversity(ttax_rpk, index="shannon")
ttax_rpk$mid<-rownames(ttax_rpk)
# read in metadata
meta_med<-read.csv(file= "./data_generated/meta_species.csv")
med_id<-subset(meta_med, select=mid)
ttax_rpk<-inner_join(med_id, ttax_rpk, by="mid")
ttax_rpkshan <- subset(ttax_rpk,select=c(mid,shannon))
meta_shan<-inner_join(meta_med, ttax_rpkshan, by="mid")

shannon.paee_pam <- lmer( shannon ~ paee_pamwk + age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+(1 | SubjectID) , data=meta_shan)
summary(shannon.paee_pam)

# categorize paee_pam
meta_shan$paee_pamwkq <- cut(meta_shan$paee_pamwk, breaks=c(-Inf, 16, 22, 30, Inf), label=c("Q1: \u226416", "Q2: 17~22", "Q3: 23~30", "Q4: >30"))
meta_shan$paee_pamwkqmean[meta_shan$paee_pamwkq=="Q1: \u226416"]<-1
meta_shan$paee_pamwkqmean[meta_shan$paee_pamwkq=="Q2: 17~22"]<-2
meta_shan$paee_pamwkqmean[meta_shan$paee_pamwkq=="Q3: 23~30"]<-3
meta_shan$paee_pamwkqmean[meta_shan$paee_pamwkq=="Q4: >30"]<-4

# boxplot for paee_pam and shannon association
png(file="./figure_generated/shannon_paee_pam_box.png",width=2000,height=1900, pointsize=70)
ggplot(meta_shan, aes(x = factor(paee_pamwkq), y = shannon, fill=meta_shan$paee_pamwkqmean)) +
  geom_boxplot(colour = "black",lwd=1, outlier.color = "black",outlier.size = 3,notch = TRUE)+theme_bw() +scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  scale_y_continuous(limits=c(0,5), breaks=seq(0,5,1))+
  xlab("Recent total PA (MET-hours/week)")+
  ylab("Shannon Diversity Index")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=2),
        plot.title = element_blank(),
        axis.title.x=element_text(color = "Black",size=70),
        axis.title.y=element_text(color = "Black",size=70),
        axis.text.y=element_text(color = "Black",size=70),
        axis.text.x=element_text(color = "Black",size=70))
dev.off()
