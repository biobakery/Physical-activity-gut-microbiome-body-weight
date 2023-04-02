
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity and abundance of A. putredinis in relation to body weight change
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
library(colourvalues)

# read in metadata
meta_med<-read.csv(file="./data_generated/meta_species.csv")
meta_med$mvpa_paqlong <- meta_med$vpa_paqlong + meta_med$mpa_paqlong

summary(meta_med$act_paqlong_sec)
summary(meta_med$act_paqlong_como)

# read in taxonomy data
tax_rpk_name <-   read.table(file = './data_generated/bugs_dna_929_unFilt.tsv',
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

tax_rpk<-tax_rpk_species+1
tax_rpk_rel <-
  sweep(tax_rpk,
        2,
        STATS = colSums(tax_rpk),
        FUN = "/")
tax_rpk_rel<-log10(tax_rpk_rel)
ttax_rel <- as.data.frame(t(tax_rpk_rel))
ttax_rel<- ttax_rel[,order(-colMeans(ttax_rel))]
ttax_rel$mid<-rownames(ttax_rel)
ttax_rel<-inner_join(meta_med, ttax_rel, by="mid")
rownames(ttax_rel)<-ttax_rel$mid

# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of tdee_pam colored by bmi_dlw
png(file="./figure_generated/pa_bmi_Alistipes.png",width=2000,height=2000, pointsize=80)
ggplot(ttax_rel, aes(s__Alistipes_putredinis.y, paee_pamwk)) + geom_point(aes(color = bmi_dlw),size=20, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  ylab("Recent total PA (MET-hours/week)")+theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.text = element_text(size=80),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=80),
        axis.text.x=element_text(color = "Black",size=80),
        axis.title.y=element_text(color = "Black",size=80),
        axis.title.x=element_blank())
dev.off()

# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of tdee_pam colored by pfat_dlw
png(file="./figure_generated/pa_pfat_Alistipes.png",width=2000,height=2000, pointsize=80)
ggplot(ttax_rel, aes(s__Alistipes_putredinis.y, paee_pamwk)) + geom_point(aes(color = pfat_dlw),size=20, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  xlab("Abundance of Alistipes putredinis")+
  ylab("Recent total PA (MET-hours/week)")+theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.text = element_text(size=80),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=80),
        axis.text.x=element_text(color = "Black",size=80),
        axis.title.y=element_text(color = "Black",size=80),
        axis.title.x=element_text(color = "Black",size=80))
dev.off()

# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of act_paqlong colored by wtchgsto21
png(file="./figure_generated/pa_sec_wtchgsto21_Alistipes.png",width=2000,height=2000, pointsize=80)
ggplot(ttax_rel, aes(s__Alistipes_putredinis.y, act_paqlong_sec)) + geom_point(aes(color = wtchgsto21),size=20, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#FDE725FF", "#B4DD2CFF","#5DC963FF","#21908CFF","#3B528BFF","#482878FF","#440154FF"))+
  ylim(0, 90)+
  xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  ylab("Long-term PA (MET-hours/week)")+theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.text = element_text(size=80),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=80),
        axis.text.x=element_text(color = "Black",size=80),
        axis.title.y=element_text(color = "Black",size=80),
        axis.title.x=element_blank())
dev.off()

