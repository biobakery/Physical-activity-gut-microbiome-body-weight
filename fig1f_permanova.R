
##############################################################################################################################################
# 1) Purposes: conduct PERMANOVA using brey-curtis to calculation the contribution from the exposure, outcome, and covariates to the variarion of overall microbial structure
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

# read in taxonomy data
tax_rpk_name <-   read.table(    file = './data_generated/bugs_dna_929_unFilt.tsv',
                                 sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)

# only keep species-level features
tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))

rownames(tax_rpk_species)<-tax_rpk_species$species
tax_rpk_species<-tax_rpk_species[,-c(1:8)]

# read in metadata
meta_med<-read.csv(file="./data_generated/meta_species.csv")
meta_med$mvpa_paqlong <- meta_med$vpa_paqlong + meta_med$mpa_paqlong
rownames(meta_med)<-meta_med$mid

med_id<-meta_med %>% select(mid) 
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

# calculate bray-curtis matrix
ttax_rpk_rel_bc <- vegdist(ttax_rpk_rel, "bray")

# PERMANOVA using adonis, strata=as.factor(meta_med$SubjectID)
x <- adonis(formula = ttax_rpk_rel_bc ~ paee_pam , strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
paee_pam <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(paee_pam)

x <- adonis(formula = ttax_rpk_rel_bc ~ mets_vigpam , strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
mets_vigpam <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(mets_vigpam)

x <- adonis(formula = ttax_rpk_rel_bc ~ mets_modpam , strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
mets_modpam <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(mets_modpam)

x <- adonis(formula = ttax_rpk_rel_bc ~ mets_ltpam , strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
mets_ltpam <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(mets_ltpam)

x <- adonis(formula = ttax_rpk_rel_bc ~ act_paqlong , strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
act_paqlong <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(act_paqlong)

x <- adonis(formula = ttax_rpk_rel_bc ~ vpa_paqlong ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
vpa_paqlong <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(vpa_paqlong)

x <- adonis(formula = ttax_rpk_rel_bc ~ mpa_paqlong ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
mpa_paqlong <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(mpa_paqlong)

x <- adonis(formula = ttax_rpk_rel_bc ~ lpa_paqlong ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
lpa_paqlong <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(lpa_paqlong)

x <- adonis(formula = ttax_rpk_rel_bc ~ bmi_dlw ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
bmi_dlw <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(bmi_dlw)

x <- adonis(formula = ttax_rpk_rel_bc ~ pfat_dlw ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
pfat_dlw <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(pfat_dlw)

x <- adonis(formula = ttax_rpk_rel_bc ~ weightchg_blood ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
weightchg_blood <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(weightchg_blood)

x <- adonis(formula = ttax_rpk_rel_bc ~ wtchgsto21 ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
wtchgsto21 <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(wtchgsto21)

x <- adonis(formula = ttax_rpk_rel_bc ~ lg2crp ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
lg2crp <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(lg2crp)

x <- adonis(formula = ttax_rpk_rel_bc ~ lg2hba1c ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
lg2hba1c <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(lg2hba1c)

x <- adonis(formula = ttax_rpk_rel_bc ~ ahei ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
ahei <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(ahei)

x <- adonis(formula = ttax_rpk_rel_bc ~  calor122cn ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
calor122cn <- x$aov.tab[2,c("R2","Pr(>F)")]
summary(calor122cn)

x <- adonis(formula = ttax_rpk_rel_bc ~  ant_12mo_qu ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
ant_12mo_qu <- x$aov.tab[2,c("R2","Pr(>F)")]
summary(ant_12mo_qu)

x <- adonis(formula = ttax_rpk_rel_bc ~ probio_2mo_qu ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
probio_2mo_qu <- x$aov.tab[3,c("R2","Pr(>F)")]
summary(probio_2mo_qu)

x <- adonis(formula = ttax_rpk_rel_bc ~ stool_type ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
stool_type <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(stool_type)

x <- adonis(formula = ttax_rpk_rel_bc ~ age_fecal ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
age_fecal <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(age_fecal)

x <- adonis(formula = ttax_rpk_rel_bc ~ smk ,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
x
smk <- x$aov.tab[1,c("R2","Pr(>F)")]
summary(smk)

all_taxonomy <-rbind(paee_pam,
                     mets_vigpam,
                     mets_modpam,
                     mets_ltpam,
                     act_paqlong, 
                     vpa_paqlong, 
                     mpa_paqlong,
                     lpa_paqlong,
                     bmi_dlw, 
                     pfat_dlw,
                     weightchg_blood,
                     wtchgsto21, 
                     lg2crp, 
                     lg2hba1c,
                     ahei,
                     calor122cn,
                     ant_12mo_qu,
                     probio_2mo_qu,
                     stool_type,
                     age_fecal,
                     smk)

component<- c('Recent total PA',
              'Recent vigorous PA',
              'Recent moderate PA',
              'Recent light PA',
              'Long-term total PA',
              'Long-term vigorous PA',
              'Long-term moderate PA',
              'Long-term light PA',
              'BMI',
              'Fat mass %',
              '6-month weight change',
              'Weight change since age 21',
              'C-reactive protein',
              'Hemoglobin A1c',
              'Alternative Healthy Eating Index',
              'Total energy intake',
              'Antibiotics use',
              'Probiotics use',
              'Bristol score',
              'Age',
              'Smoking')  

all_taxonomy<-cbind(all_taxonomy, component)
all_taxonomy_pa<-all_taxonomy[1:8,]
all_taxonomy_pa$cat<-"PA"

all_taxonomy_wt<-all_taxonomy[9:12,]
all_taxonomy_wt$cat<-"Body weight"

all_taxonomy_bio<-all_taxonomy[13:14,]
all_taxonomy_bio$cat<-"Plasma biomarkers"

all_taxonomy_cov<-all_taxonomy[15:21,]
all_taxonomy_cov$cat<-"Covariables"

all_taxonomy_tax<-rbind(all_taxonomy_pa, all_taxonomy_wt, all_taxonomy_bio, all_taxonomy_cov)
all_taxonomy_tax$cat<-as.factor(all_taxonomy_tax$cat)

color_code<-c("PA"="#3EAB07",
              "Body weight"="#FF1713",
              "Plasma biomarkers"="#FFB713",
              "Covariables"="#8dd3c7")

all_taxonomy_tax$stars <- cut(all_taxonomy_tax$`Pr(>F)`, breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))
write.csv(all_taxonomy_tax, file="./data_generated/permanova_tax.csv")

all_taxonomy_tax<-read.csv(file="./data_generated/permanova_tax.csv")

level_y_order <- factor(all_taxonomy_tax$component, level = rev(all_taxonomy_tax$component))

# create figure to show PERMANOVA results
png(file="./figure_generated/fig2b_barplot_permanova.png",width=2600, height=2000, pointsize=50)
ggplot(all_taxonomy_tax, aes(level_y_order, R2, fill=cat)) +
  geom_bar(stat="identity")+
  theme_light()+scale_y_continuous(labels=percent, limits = c(0, 0.01))+
  geom_text(aes(label=stars), color="black", size=40) +
  scale_fill_manual(values = color_code)+ coord_flip()+
  theme(axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.position = "none",
        legend.text = element_text(size = 60,color="black"),
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=60,color="black"),
        axis.text.x = element_text(size=60,color="black"))  
dev.off()
