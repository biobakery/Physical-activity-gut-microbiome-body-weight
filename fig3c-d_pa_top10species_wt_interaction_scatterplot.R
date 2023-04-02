
##############################################################################################################################################
# 1) Purposes: creating scatterplot of interaction analysis between physical activity and the top 10 abundant species in relation to body weight change
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
meta_species$s__Eubacterium_rectale2c <- with(meta_species, ifelse(s__Eubacterium_rectale> 6.616, 1, 0))

meta_species$cutx[meta_species$s__Eubacterium_rectale2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_rectale2c==1]<-"> median"
png(file="./figure_generated/e_rectale_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium rectale")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#2. s__Faecalibacterium_prausnitzii
summary(meta_species$s__Faecalibacterium_prausnitzii)
meta_species$s__Faecalibacterium_prausnitzii2c <- with(meta_species, ifelse(s__Faecalibacterium_prausnitzii> 6.452, 1, 0))

meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==1]<-"> median"
png(file="./figure_generated/f_prausnitzii_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Faecalibacterium prausnitzii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#3. s__Subdoligranulum_unclassified
summary(meta_species$s__Subdoligranulum_unclassified)
meta_species$s__Subdoligranulum_unclassified2c <- with(meta_species, ifelse(s__Subdoligranulum_unclassified> 5.234, 1, 0))

meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==1]<-"> median"
png(file="./figure_generated/s_unclassified_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Subdoligranulum unclassified")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#4. s__Bacteroides_uniformis
summary(meta_species$s__Bacteroides_uniformis)
meta_species$s__Bacteroides_uniformis2c <- with(meta_species, ifelse(s__Bacteroides_uniformis> 4.059, 1, 0))

meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==1]<-"> median"
png(file="./figure_generated/b_uniformis_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides uniformis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#5. s__Prevotella_copri
summary(meta_species$s__Prevotella_copri)
meta_species$s__Prevotella_copri2c <- with(meta_species, ifelse(s__Prevotella_copri> 0, 1, 0))

meta_species$cutx[meta_species$s__Prevotella_copri2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Prevotella_copri2c==1]<-"> median"
png(file="./figure_generated/p_copri_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Prevotella copri")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#6. s__Alistipes_putredinis
summary(meta_species$s__Alistipes_putredinis)
meta_species$s__Alistipes_putredinis2c <- with(meta_species, ifelse(s__Alistipes_putredinis> 3.0033, 1, 0))

meta_species$cutx[meta_species$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/a_putredinis_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#7. s__Eubacterium_siraeum
summary(meta_species$s__Eubacterium_siraeum)
meta_species$s__Eubacterium_siraeum2c <- with(meta_species, ifelse(s__Eubacterium_siraeum> 0.2456, 1, 0))

meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==1]<-"> median"
png(file="./figure_generated/e_siraeum_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium siraeum")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#8. s__Ruminococcus_bromii
summary(meta_species$s__Ruminococcus_bromii)
meta_species$s__Ruminococcus_bromii2c <- with(meta_species, ifelse(s__Ruminococcus_bromii> 0.6912, 1, 0))

meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==1]<-"> median"
png(file="./figure_generated/r_bromii_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Ruminococcus bromii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#9. s__Bacteroides_vulgatus
summary(meta_species$s__Bacteroides_vulgatus)
meta_species$s__Bacteroides_vulgatus2c <- with(meta_species, ifelse(s__Bacteroides_vulgatus> 1.5450, 1, 0))

meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==1]<-"> median"
png(file="./figure_generated/b_vulgatus_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides vulgatus")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#10. s__Bacteroides_stercoris
summary(meta_species$s__Bacteroides_stercoris)
meta_species$s__Bacteroides_stercoris2c <- with(meta_species, ifelse(s__Bacteroides_stercoris> 0.02914, 1, 0))

meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==1]<-"> median"
png(file="./figure_generated/b_stercoris_paee_pam_bmi_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, bmi_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides stercoris")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~(kg/(m^2))))+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

####### pfat_dlw ########

#1. s__Eubacterium_rectale

meta_species$cutx[meta_species$s__Eubacterium_rectale2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_rectale2c==1]<-"> median"
png(file="./figure_generated/e_rectale_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium rectale")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#2. s__Faecalibacterium_prausnitzii

meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==1]<-"> median"
png(file="./figure_generated/f_prausnitzii_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Faecalibacterium prausnitzii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#3. s__Subdoligranulum_unclassified

meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==1]<-"> median"
png(file="./figure_generated/s_unclassified_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Subdoligranulum unclassified")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#4. s__Bacteroides_uniformis

meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==1]<-"> median"
png(file="./figure_generated/b_uniformis_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides uniformis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#5. s__Prevotella_copri

meta_species$cutx[meta_species$s__Prevotella_copri2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Prevotella_copri2c==1]<-"> median"
png(file="./figure_generated/p_copri_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Prevotella copri")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#6. s__Alistipes_putredinis

meta_species$cutx[meta_species$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/a_putredinis_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#7. s__Eubacterium_siraeum

meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==1]<-"> median"
png(file="./figure_generated/e_siraeum_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium siraeum")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#8. s__Ruminococcus_bromii

meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==1]<-"> median"
png(file="./figure_generated/r_bromii_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Ruminococcus bromii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#9. s__Bacteroides_vulgatus

meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==1]<-"> median"
png(file="./figure_generated/b_vulgatus_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides vulgatus")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#10. s__Bacteroides_stercoris

meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==1]<-"> median"
png(file="./figure_generated/b_stercoris_paee_pam_pfat_dlw.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, pfat_dlw, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides stercoris")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("pfat_dlw")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

####### wtchgsto21 ########

#1. s__Eubacterium_rectale

meta_species$cutx[meta_species$s__Eubacterium_rectale2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_rectale2c==1]<-"> median"
png(file="./figure_generated/e_rectale_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
 # scale_y_continuous(breaks=seq(-30,50,10))+
#  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium rectale")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#2. s__Faecalibacterium_prausnitzii

meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==1]<-"> median"
png(file="./figure_generated/f_prausnitzii_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  # scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Faecalibacterium prausnitzii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#3. s__Subdoligranulum_unclassified

meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==1]<-"> median"
png(file="./figure_generated/s_unclassified_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  # scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Subdoligranulum unclassified")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#4. s__Bacteroides_uniformis

meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==1]<-"> median"
png(file="./figure_generated/b_uniformis_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  # scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides uniformis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#5. s__Prevotella_copri

meta_species$cutx[meta_species$s__Prevotella_copri2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Prevotella_copri2c==1]<-"> median"
png(file="./figure_generated/p_copri_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  # scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Prevotella copri")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#6. s__Alistipes_putredinis

meta_species$cutx[meta_species$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/a_putredinis_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  # scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#7. s__Eubacterium_siraeum

meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==1]<-"> median"
png(file="./figure_generated/e_siraeum_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  #  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium siraeum")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#8. s__Ruminococcus_bromii

meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==1]<-"> median"
png(file="./figure_generated/r_bromii_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  #  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Ruminococcus bromii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#9. s__Bacteroides_vulgatus

meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==1]<-"> median"
png(file="./figure_generated/b_vulgatus_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  # scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides vulgatus")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#10. s__Bacteroides_stercoris

meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==1]<-"> median"
png(file="./figure_generated/b_stercoris_act_paqlong_wtchgsto21.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act_paqlong, wtchgsto21, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  # scale_y_continuous(breaks=seq(-30,50,10))+
  # scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides stercoris")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("wtchgsto21")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

####### weightchg_blood ########

#1. s__Eubacterium_rectale

summary(meta_species_12stool$s__Eubacterium_rectale)
table(meta_species_12stool$s__Eubacterium_rectale)

graphics.off()
histogram(meta_species_12stool$s__Eubacterium_rectale)

meta_species_12stool$s__Eubacterium_rectale2c <- with(meta_species_12stool, ifelse(s__Eubacterium_rectale> 6.518, 1, 0))
table(meta_species_12stool$s__Eubacterium_rectale2c)
#0   1 
#209 208

meta_species_12stool$cutx[meta_species_12stool$s__Eubacterium_rectale2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Eubacterium_rectale2c==1]<-"> median"
png(file="./figure_generated/e_rectale_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium rectale")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#2. s__Faecalibacterium_prausnitzii
summary(meta_species_12stool$s__Faecalibacterium_prausnitzii)
table(meta_species_12stool$s__Faecalibacterium_prausnitzii)

graphics.off()
histogram(meta_species_12stool$s__Faecalibacterium_prausnitzii)

meta_species_12stool$s__Faecalibacterium_prausnitzii2c <- with(meta_species_12stool, ifelse(s__Faecalibacterium_prausnitzii> 6.452, 1, 0))
table(meta_species_12stool$s__Faecalibacterium_prausnitzii2c)
#0   1 
#209 208
meta_species_12stool$cutx[meta_species_12stool$s__Faecalibacterium_prausnitzii2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Faecalibacterium_prausnitzii2c==1]<-"> median"
png(file="./figure_generated/f_prausnitzii_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Faecalibacterium prausnitzii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#3. s__Subdoligranulum_unclassified
summary(meta_species_12stool$s__Subdoligranulum_unclassified)
table(meta_species_12stool$s__Subdoligranulum_unclassified)

graphics.off()
histogram(meta_species_12stool$s__Subdoligranulum_unclassified)

meta_species_12stool$s__Subdoligranulum_unclassified2c <- with(meta_species_12stool, ifelse(s__Subdoligranulum_unclassified> 5.037, 1, 0))
table(meta_species_12stool$s__Subdoligranulum_unclassified2c)
#0   1 
#208 209
meta_species_12stool$cutx[meta_species_12stool$s__Subdoligranulum_unclassified2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Subdoligranulum_unclassified2c==1]<-"> median"
png(file="./figure_generated/s_unclassified_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Subdoligranulum unclassified")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#4. s__Bacteroides_uniformis
summary(meta_species_12stool$s__Bacteroides_uniformis)
table(meta_species_12stool$s__Bacteroides_uniformis)

graphics.off()
histogram(meta_species_12stool$s__Bacteroides_uniformis)

meta_species_12stool$s__Bacteroides_uniformis2c <- with(meta_species_12stool, ifelse(s__Bacteroides_uniformis> 4.059, 1, 0))
table(meta_species_12stool$s__Bacteroides_uniformis2c)
#0   1 
#208 209
meta_species_12stool$cutx[meta_species_12stool$s__Bacteroides_uniformis2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Bacteroides_uniformis2c==1]<-"> median"
png(file="./figure_generated/b_uniformis_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides uniformis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#5. s__Prevotella_copri
summary(meta_species_12stool$s__Prevotella_copri)
table(meta_species_12stool$s__Prevotella_copri)

graphics.off()
histogram(meta_species_12stool$s__Prevotella_copri)

meta_species_12stool$s__Prevotella_copri2c <- with(meta_species_12stool, ifelse(s__Prevotella_copri> 0, 1, 0))
table(meta_species_12stool$s__Prevotella_copri2c)
#0   1 
#326 91
meta_species_12stool$cutx[meta_species_12stool$s__Prevotella_copri2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Prevotella_copri2c==1]<-"> median"
png(file="./figure_generated/p_copri_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Prevotella copri")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#6. s__Alistipes_putredinis
summary(meta_species_12stool$s__Alistipes_putredinis)
table(meta_species_12stool$s__Alistipes_putredinis)

graphics.off()
histogram(meta_species_12stool$s__Alistipes_putredinis)

meta_species_12stool$s__Alistipes_putredinis2c <- with(meta_species_12stool, ifelse(s__Alistipes_putredinis> 3.1016, 1, 0))
table(meta_species_12stool$s__Alistipes_putredinis2c)
#0   1 
#209 208
meta_species_12stool$cutx[meta_species_12stool$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/a_putredinis_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#7. s__Eubacterium_siraeum
summary(meta_species_12stool$s__Eubacterium_siraeum)
table(meta_species_12stool$s__Eubacterium_siraeum)

graphics.off()
histogram(meta_species_12stool$s__Eubacterium_siraeum)

meta_species_12stool$s__Eubacterium_siraeum2c <- with(meta_species_12stool, ifelse(s__Eubacterium_siraeum> 0.2039, 1, 0))
table(meta_species_12stool$s__Eubacterium_siraeum2c)
#0   1 
#208 209
meta_species_12stool$cutx[meta_species_12stool$s__Eubacterium_siraeum2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Eubacterium_siraeum2c==1]<-"> median"
png(file="./figure_generated/e_siraeum_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium siraeum")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#8. s__Ruminococcus_bromii
summary(meta_species_12stool$s__Ruminococcus_bromii)
table(meta_species_12stool$s__Ruminococcus_bromii)

graphics.off()
histogram(meta_species_12stool$s__Ruminococcus_bromii)

meta_species_12stool$s__Ruminococcus_bromii2c <- with(meta_species_12stool, ifelse(s__Ruminococcus_bromii> 0.7604, 1, 0))
table(meta_species_12stool$s__Ruminococcus_bromii2c)
#0   1 
#209 208
meta_species_12stool$cutx[meta_species_12stool$s__Ruminococcus_bromii2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Ruminococcus_bromii2c==1]<-"> median"
png(file="./figure_generated/r_bromii_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Ruminococcus bromii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#9. s__Bacteroides_vulgatus
summary(meta_species_12stool$s__Bacteroides_vulgatus)
table(meta_species_12stool$s__Bacteroides_vulgatus)

graphics.off()
histogram(meta_species_12stool$s__Bacteroides_vulgatus)

meta_species_12stool$s__Bacteroides_vulgatus2c <- with(meta_species_12stool, ifelse(s__Bacteroides_vulgatus> 1.544, 1, 0))
table(meta_species_12stool$s__Bacteroides_vulgatus2c)
#0   1 
#208 209
meta_species_12stool$cutx[meta_species_12stool$s__Bacteroides_vulgatus2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Bacteroides_vulgatus2c==1]<-"> median"
png(file="./figure_generated/b_vulgatus_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides vulgatus")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#10. s__Bacteroides_stercoris
summary(meta_species_12stool$s__Bacteroides_stercoris)
table(meta_species_12stool$s__Bacteroides_stercoris)

graphics.off()
histogram(meta_species_12stool$s__Bacteroides_stercoris)

meta_species_12stool$s__Bacteroides_stercoris2c <- with(meta_species_12stool, ifelse(s__Bacteroides_stercoris> 0.00672, 1, 0))
table(meta_species_12stool$s__Bacteroides_stercoris2c)
#0   1 
#209 208
meta_species_12stool$cutx[meta_species_12stool$s__Bacteroides_stercoris2c==0]<-"\u2264 median"
meta_species_12stool$cutx[meta_species_12stool$s__Bacteroides_stercoris2c==1]<-"> median"
png(file="./figure_generated/b_stercoris_paee_pam_weightchg_blood.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species_12stool, aes(paee_pamwk, weightchg_blood, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  scale_y_continuous(breaks=seq(-10,10,2.5))+
  scale_x_continuous(breaks=seq(0,160,20))+
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides stercoris")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("weightchg_blood")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

####### lg2hba1c ########

#1. s__Eubacterium_rectale

meta_species$cutx[meta_species$s__Eubacterium_rectale2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_rectale2c==1]<-"> median"
png(file="./figure_generated/e_rectale_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium rectale")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#2. s__Faecalibacterium_prausnitzii

meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==1]<-"> median"
png(file="./figure_generated/f_prausnitzii_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Faecalibacterium prausnitzii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#3. s__Subdoligranulum_unclassified

meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==1]<-"> median"
png(file="./figure_generated/s_unclassified_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Subdoligranulum unclassified")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#4. s__Bacteroides_uniformis

meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==1]<-"> median"
png(file="./figure_generated/b_uniformis_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides uniformis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#5. s__Prevotella_copri

meta_species$cutx[meta_species$s__Prevotella_copri2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Prevotella_copri2c==1]<-"> median"
png(file="./figure_generated/p_copri_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Prevotella copri")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#6. s__Alistipes_putredinis

meta_species$cutx[meta_species$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/a_putredinis_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#7. s__Eubacterium_siraeum

meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==1]<-"> median"
png(file="./figure_generated/e_siraeum_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium siraeum")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#8. s__Ruminococcus_bromii

meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==1]<-"> median"
png(file="./figure_generated/r_bromii_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Ruminococcus bromii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#9. s__Bacteroides_vulgatus

meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==1]<-"> median"
png(file="./figure_generated/b_vulgatus_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides vulgatus")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#10. s__Bacteroides_stercoris

meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==1]<-"> median"
png(file="./figure_generated/b_stercoris_paee_pam_lg2hba1c.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2hba1c, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides stercoris")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2hba1c")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

####### lg2crp ########

#1. s__Eubacterium_rectale

meta_species$cutx[meta_species$s__Eubacterium_rectale2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_rectale2c==1]<-"> median"
png(file="./figure_generated/e_rectale_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium rectale")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#2. s__Faecalibacterium_prausnitzii

meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Faecalibacterium_prausnitzii2c==1]<-"> median"
png(file="./figure_generated/f_prausnitzii_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Faecalibacterium prausnitzii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#3. s__Subdoligranulum_unclassified

meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Subdoligranulum_unclassified2c==1]<-"> median"
png(file="./figure_generated/s_unclassified_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Subdoligranulum unclassified")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#4. s__Bacteroides_uniformis

meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_uniformis2c==1]<-"> median"
png(file="./figure_generated/b_uniformis_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides uniformis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#5. s__Prevotella_copri

meta_species$cutx[meta_species$s__Prevotella_copri2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Prevotella_copri2c==1]<-"> median"
png(file="./figure_generated/p_copri_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Prevotella copri")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#6. s__Alistipes_putredinis

meta_species$cutx[meta_species$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/a_putredinis_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#7. s__Eubacterium_siraeum

meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Eubacterium_siraeum2c==1]<-"> median"
png(file="./figure_generated/e_siraeum_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Eubacterium siraeum")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#8. s__Ruminococcus_bromii

meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Ruminococcus_bromii2c==1]<-"> median"
png(file="./figure_generated/r_bromii_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Ruminococcus bromii")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#9. s__Bacteroides_vulgatus

meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_vulgatus2c==1]<-"> median"
png(file="./figure_generated/b_vulgatus_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides vulgatus")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

#10. s__Bacteroides_stercoris

meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Bacteroides_stercoris2c==1]<-"> median"
png(file="./figure_generated/b_stercoris_paee_pam_lg2crp.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(paee_pamwk, lg2crp, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Bacteroides stercoris")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("lg2crp")+theme_classic()+
  theme(title = element_text(size = 50),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()
