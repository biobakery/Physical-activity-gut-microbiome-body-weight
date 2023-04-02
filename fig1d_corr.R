##############################################################################################################################################
# 1) Purposes: conduct correlation analysis among the main exposure and outcome variables, 
#              including short- and long-term physical activity and BMI, fat mass%, short- and long-term weight change, and plasma HbA1c and CRP
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
library("Hmisc")
library(reshape2)
library(expss)

meta_species<-read.csv(file="./data_generated/meta_species.csv", header = TRUE)

cor_data <- subset(meta_species, select=c(act_paqlong,bmi_dlw,wtchgsto21,calorie,alco_a,prot_a,aprot_a,vprot_a,dprot_a,fat_a,afat_a,vfat_a,dfat_a,carbo_a,aofib_a))
cor <- rcorr(as.matrix(cor_data ),type=c("spearman"))
cor

meta_species$'Recent_total_PA' <- meta_species$paee_pamwk
meta_species$'Long_term_total_PA' <- meta_species$act_paqlong
meta_species$'BMI' <- meta_species$bmi_dlw
meta_species$'Fat_mass_percentage' <- meta_species$pfat_dlw
meta_species$'Six_month_weight_change' <- meta_species$bodyweightchg_dlw
meta_species$'Weight_change_since_age_21' <- meta_species$wtchgsto21
meta_species$'Plasma_CRP' <- meta_species$lgcrp
meta_species$'Plasma_HbA1c' <- meta_species$lghba1c

pawt_data <- subset(meta_species, select=c('Recent_total_PA', 
                                           'Long_term_total_PA', 
                                           'BMI', 
                                           'Fat_mass_percentage', 
                                           'Six_month_weight_change', 
                                           'Weight_change_since_age_21', 
                                           'Plasma_CRP', 
                                           'Plasma_HbA1c'))

cormat <- rcorr(as.matrix(pawt_data),type=c("spearman"))
cormat
rcx <- cormat
str(rcx)
df.rcx<-data.frame(rcx$r)
df.rcx<-round(df.rcx,2)
df.rcx

# Get lower triangle of the correlation matrix
get_lower_tri<-function(df.rcx){
  df.rcx[lower.tri(df.rcx)] <- NA
  return(df.rcx)
}
lower_tri <- get_lower_tri(df.rcx)
lower_tri
lower_tri$var1<-rownames(lower_tri)

melted_cormat <- melt(lower_tri, id.vars='var1', na.rm = TRUE)
melted_cormat
# Heatmap

level_x_order <- factor(melted_cormat$variable, level = c('Recent_total_PA', 
                                                          'Long_term_total_PA', 
                                                          'BMI', 
                                                          'Fat_mass_percentage', 
                                                          'Six_month_weight_change', 
                                                          'Weight_change_since_age_21', 
                                                          'Plasma_CRP', 
                                                          'Plasma_HbA1c'))

level_y_order <- factor(melted_cormat$var1, level = c('Recent_total_PA', 
                                                      'Long_term_total_PA', 
                                                      'BMI', 
                                                      'Fat_mass_percentage', 
                                                      'Six_month_weight_change', 
                                                      'Weight_change_since_age_21', 
                                                      'Plasma_CRP', 
                                                      'Plasma_HbA1c'))

ggheatmap <- ggplot(melted_cormat, aes(level_x_order, level_y_order, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  geom_text(aes(label=value), color="black", size=30, show.legend = TRUE) +
  theme_minimal()+ # minimal theme
  theme(legend.title = element_text(size = 80,color="black"),
        legend.text = element_text(size = 80,color="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 80, hjust = 1,color="black"),
        axis.text.y = element_text(size = 80,color="black"),
        axis.title.x= element_blank(),
        axis.title.y= element_blank())+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

png(file="./figure_generated/fig1_correlation_pawt.png",width=3500,height=3500, pointsize=50)
grid.arrange(ggheatmap, ncol=1, nrow=1)
dev.off()


