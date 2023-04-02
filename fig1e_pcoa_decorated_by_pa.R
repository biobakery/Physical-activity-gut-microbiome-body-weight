##############################################################################################################################################
# 1) Purposes: conduct PCo analysis on species level taxonomy data and plot PCo1 and 2 colored by physical activity
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

library(tidyverse)
library(vegan)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)

# read in taxonomy data
tax_rpk_name <-   read.table( "/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/taxonomy/bugs_dna_929_unFilt.tsv",
                                 sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)

# only keep species-level features
tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))
rownames(tax_rpk_species)<-tax_rpk_species$species
tax_rpk_species<-tax_rpk_species[,-c(1:8)]

# read in metadata
meta_med<-read.csv(file= "./data_generated/meta_ffq_8612_biomarker.csv")

med_id<-subset(meta_med, select=mid)

ttax_rpk <- as.data.frame(t(tax_rpk_species))

ttax_rpk$mid<-rownames(ttax_rpk)
ttax_rpk<-inner_join(med_id, ttax_rpk, by="mid")

rownames(ttax_rpk)<-ttax_rpk$mid
ttax_rpk <- subset(ttax_rpk, select = -mid)

ttax_rpk_rel <-
  sweep(ttax_rpk,
        1,
        STATS = rowSums(ttax_rpk),
        FUN = "/")

# calculate bray-crutis matrix
ttax_rpk_rel_bc <- vegdist(ttax_rpk_rel, "bray")

# PCOA using cmdscale
mod <- cmdscale(ttax_rpk_rel_bc, eig = T)
pcoap <- data.frame(mod$points)
pcoap$mid<-rownames(pcoap)

# percentages of variation explained by PCO1 & 2
mod$eig[1]/sum(mod$eig)
mod$eig[2]/sum(mod$eig)

med_score<-subset(meta_med, select = c(mid, act_paqq))
pcoap<-inner_join(pcoap, med_score, by="mid")
ttax_rpk_rel$mid<-rownames(ttax_rpk_rel)
ttax_pcoap<-inner_join(pcoap, ttax_rpk_rel, by="mid")
write.csv(ttax_pcoap, file="./pcoa_species2.csv")
write.csv(pcoap, file="./pcoa_species.csv")

# PCOA figure colored by act_paqq
png(
  file = "./pcoa_act_paqq.png",
  width = 2400,
  height = 2400,
  pointsize = 50
)
ggplot(data = pcoap, aes(X1, X2, colour = pcoap$act_paqq
                         )) +
  geom_point(size=22, alpha=0.8) + 
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +
  theme_classic()+
  theme(legend.position = "bottom",
    legend.title = element_text(size = 80),
    legend.text = element_text(size = 80),
    axis.title.x = element_text(size = 80),
    axis.title.y = element_text(size = 80),
    axis.text.y = element_text(size = 80),
    axis.text.x = element_text(size = 80),
    axis.line = element_line(color = "Black", size=3)
  ) +
scale_color_viridis(name = "Physical activity Met-hours/week")
dev.off()

meta_species<-read.csv(file="./data_generated/meta_species.csv", header = TRUE)
summary(meta_species$paee_pamwk)

# PCOA figure colored by physical activity: act_paqlong
png(file = "./figure_generated/pcoa_by_act_paqlong.png",width = 2400,height = 2400,pointsize = 50)
ggplot(data = meta_species, aes(X1, X2)) + geom_point(aes(color = act_paqlong),size=22, alpha=0.8, position = position_jitter()) +
  scale_color_gradientn(colours = c("#FDE725FF", "#B4DD2CFF","#5DC963FF","#21908CFF","#3B528BFF","#482878FF","#440154FF"))+
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 80),
        legend.text = element_text(size = 70),
        axis.title.x = element_text(size = 80),
        axis.title.y = element_text(size = 80),
        axis.text.y = element_text(size = 80),
        axis.text.x = element_text(size = 80),
        axis.line = element_line(color = "Black", size=3)
  )
dev.off()

# PCOA figure colored by physical activity: paee_pam
png(file = "./figure_generated/pcoa_by_paee_pamwk.png",width = 2400,height = 2400,pointsize = 50)
ggplot(data = meta_species, aes(X1, X2)) + geom_point(aes(color = paee_pamwk),size=22, alpha=0.8, position = position_jitter()) +
  scale_color_gradientn(colours = c("#FDE725FF", "#B4DD2CFF","#5DC963FF","#21908CFF","#3B528BFF","#482878FF","#440154FF"))+
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 80),
        legend.text = element_text(size = 70),
        axis.title.x = element_text(size = 80),
        axis.title.y = element_text(size = 80),
        axis.text.y = element_text(size = 80),
        axis.text.x = element_text(size = 80),
        axis.line = element_line(color = "Black", size=3)
  )
dev.off()
