##############################################################################################################################################
# 1) Purposes: creat boxplot to show the physica activity quintile in relation to body weight change in those with high and low abundance of A. putredinis
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

meta_species<-read.csv(file="./data_generated/meta_species.csv", header = TRUE)

# categorize bmi and fat percent for creating total pa against bmi and fat percent boxplot
meta_species$act_paqlongq5 <- cut(meta_species$act_paqlong_sec, breaks=c(-Inf, 19.72, 29.62, 38.43, 56.66, Inf), label=c("\u226418", "19~28", "29~37", "38~56",">56"))
meta_species$act_paqlongq5mean[meta_species$act_paqlongq5=="\u226418"]<-1
meta_species$act_paqlongq5mean[meta_species$act_paqlongq5=="19~28"]<-2
meta_species$act_paqlongq5mean[meta_species$act_paqlongq5=="29~37"]<-3
meta_species$act_paqlongq5mean[meta_species$act_paqlongq5=="38~56"]<-4
meta_species$act_paqlongq5mean[meta_species$act_paqlongq5==">56"]<-5

meta_species$s__Alistipes_putredinis2c <- with(meta_species, ifelse(s__Alistipes_putredinis> 3.0033, 1, 0))
table(meta_species$s__Alistipes_putredinis2c)

#### wtchgsto21: x-act_paqlong, y-wtchgsto21 ####
# boxplot for act_paqlong and wtchgsto21 association
png(file="./figure_generated/box_act_paqlong_wtchgsto21.png",width=500,height=1900, pointsize=70)
ggplot(meta_species, aes(x = factor(act_paqlongq5), y = wtchgsto21, fill=act_paqlongq5mean)) +
  geom_boxplot(colour = "black",lwd=1, outlier.color = "black",outlier.size = 3)+theme_bw() +scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("Weight change since age 21 (kg)")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=2),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=70),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color = "Black",size=70),
        axis.text.x=element_text(color = "Black",size=70,angle = 300,hjust = 0,vjust=1))
dev.off()
# boxplot for act_paqlong and wtchgsto21 association in s__Alistipes_putredinis2c below
png(file="./figure_generated/box_act_paqlong_wtchgsto21_lowAP2.png",width=650,height=1900, pointsize=70)
ggplot(meta_species[meta_species$s__Alistipes_putredinis2c==0, ], aes(x = factor(act_paqlongq5), y = wtchgsto21, fill=act_paqlongq5mean)) +
  geom_boxplot(colour = "black",lwd=1, outlier.color = "black",outlier.size = 3)+theme_bw() +scale_fill_gradientn(colours = c("#F5A9A9","#F78181","#FA5858","#FF0000", "#B40404"))+
  ylab("Weight change from age 21 to stool collection (kg)")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=2),
        plot.title = element_blank(),
        axis.title.y=element_text(color = "Black",size=70),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color = "Black",size=70),
        axis.text.x=element_text(color = "Black",size=70,angle = 300,hjust = 0,vjust=1))
dev.off()
# boxplot for act_paqlong and wtchgsto21 association in s__Alistipes_putredinis2c high
png(file="./figure_generated/box_act_paqlong_wtchgsto21_highAP2.png",width=500,height=1900, pointsize=70)
ggplot(meta_species[meta_species$s__Alistipes_putredinis2c==1, ], aes(x = factor(act_paqlongq5), y = wtchgsto21, fill=act_paqlongq5mean)) +
  geom_boxplot(colour = "black",lwd=1, outlier.color = "black",outlier.size = 3)+theme_bw() +scale_fill_gradientn(colours = c("#A9F5A9","#00FF00","#01DF01","#04B404", "#0B610B"))+
  ylab("Weight change since age 21 (kg)")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=2),
        plot.title = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(color = "Black",size=70,angle = 300,hjust = 0,vjust=1))
dev.off()