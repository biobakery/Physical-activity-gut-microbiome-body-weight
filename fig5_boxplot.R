
##############################################################################################################################################
# 1) Purposes: creating boxplot with recent physical activity as x-axis and body weight change from 2013-2017 as y-axis in those with high and low Allistipes putredinis
# 2) Study design:  Cross-sectional design
# 3) Endpoints: primary: Body weight change from 2013 to 2017
#               secondary: BMI at 2013 and body weight change from age 18 to stool collection
# 4) Exposures: primary: physical activity level at 2013
#               secondary: cumulative average physical activity level from 1991 to 2013
# 5) Covariates: age, diet quality, total energy intake, smoking status, antibiotic use, probiotic use, stool type
# 6) Follow-up: 2013 (Mind-Body Study)  
#############################################################################################################################################

rm(list=ls())
getwd()
setwd("/udd/nhkwa/mbs")
options(max.print=1000000)

library(viridis)
library(ggplot2)
library(grid)
library(tidyverse)
library(cowplot)
library(car)

meta_species<-read.csv(file="./data_generated/mbs_meta_species.csv", header = TRUE)

meta_species$act13mq4 <- cut(meta_species$act13m, breaks=c(-Inf, 15.2, 32.6, 54.9, Inf), label=c("\u226415", "16~33", "34~55",">55"))
meta_species$act13mq5mean[meta_species$act13mq4=="\u226415"]<-1
meta_species$act13mq5mean[meta_species$act13mq4=="16~33"]<-2
meta_species$act13mq5mean[meta_species$act13mq4=="34~55"]<-3
meta_species$act13mq5mean[meta_species$act13mq4==">55"]<-4

#s__Alistipes_putredinis
summary(meta_species$s__Alistipes_putredinis)
meta_species$s__Alistipes_putredinis2c <- with(meta_species, ifelse(s__Alistipes_putredinis> 4.5867, 1, 0))
table(meta_species$s__Alistipes_putredinis2c)

# boxplot for act13m(5g) and wtchg1317 association in s__Alistipes_putredinis2c below and high
png(file="./figure_generated/box_act13m4g_wtchg1317_lowAP.png",width=500,height=1500, pointsize=70)
ggplot(meta_species3[meta_species3$s__Alistipes_putredinis2c==1, ], aes(x = factor(act13mq4), y = wtchg1317, fill=act13mq4mean)) +
  geom_boxplot(colour = "black",lwd=1, outlier.color = "black",outlier.size = 3)+theme_bw() +scale_fill_gradientn(colours = c("#F5A9A9","#F78181","#FA5858","#FF0000", "#B40404"))+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=2),
        plot.title = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color = "Black",size=70),
        axis.text.x=element_text(color = "Black",size=70,angle = 300,hjust = 0,vjust=1))
dev.off()
png(file="./figure_generated/box_act13m4g_wtchg1317_highAP.png",width=430,height=1500, pointsize=70)
ggplot(meta_species3[meta_species3$s__Alistipes_putredinis2c==0, ], aes(x = factor(act13mq4), y = wtchg1317, fill=act13mq4mean)) +
  geom_boxplot(colour = "black",lwd=1, outlier.color = "black",outlier.size = 3)+theme_bw() +scale_fill_gradientn(colours = c("#A9F5A9","#00FF00","#01DF01","#04B404", "#0B610B"))+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=2),
        plot.title = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(color = "Black",size=70,angle = 300,hjust = 0,vjust=1))
dev.off()

