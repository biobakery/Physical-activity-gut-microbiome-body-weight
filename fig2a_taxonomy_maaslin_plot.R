
##############################################################################################################################################
# 1) Purposes: creat heatmap of the assocaitons of measures of physical activity and adiposity with abundances of per microbial species from MaAsLin regression, and prepare data for phylogenetic tree creation in GraPhlan
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
meta_med$paee_pamwk <- meta_med$paee_pam*7

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
meta_med$totMETs_paq<-as.numeric(as.character(meta_med$totMETs_paq))

# combine results form maasline regressions across different dietary variables
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}

sigpaee_pam <- sig_results(dir="./data_generated/maaslin_results/paee_pam_taxonomy/significant_results.tsv",  exposure = "paee_pam")
sigact_paqlong <- sig_results(dir="./data_generated/maaslin_results/act_paqlong_taxonomy/significant_results.tsv",  exposure = "act_paqlong")
sigbmi_dlw <- sig_results(dir="./data_generated/maaslin_results/bmi_dlw_taxonomy/significant_results.tsv",  exposure = "bmi_dlw")
sigpfat_dlw <- sig_results(dir="./data_generated/maaslin_results/pfat_dlw_taxonomy/significant_results.tsv",  exposure = "pfat_dlw")
sigweightchg_blood <- sig_results(dir="./data_generated/maaslin_results/weightchg_blood_taxonomy/significant_results.tsv",  exposure = "weightchg_blood")
sigwtchgsto21 <- sig_results(dir="./data_generated/maaslin_results/wtchgsto21_taxonomy/significant_results.tsv",  exposure = "wtchgsto21")

all_results <- function(dir, exposure, label){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}

sigpaee_pamall <- all_results(dir="./data_generated/maaslin_results/paee_pam_taxonomy/all_results.tsv",  exposure = "paee_pam", label="Recent total PA")
sigact_paqlongall <- all_results(dir="./data_generated/maaslin_results/act_paqlong_taxonomy/all_results.tsv",     exposure = "act_paqlong", label="Long-term total PA")
sigbmi_dlwall <- all_results(dir="./data_generated/maaslin_results/bmi_dlw_taxonomy/all_results.tsv",  exposure = "bmi_dlw", label="BMI")
sigpfat_dlwall <- all_results(dir="./data_generated/maaslin_results/pfat_dlw_taxonomy/all_results.tsv",    exposure = "pfat_dlw", label="Fat mass %")
sigweightchg_bloodall <- all_results(dir="./data_generated/maaslin_results/weightchg_blood_taxonomy/all_results.tsv",    exposure = "weightchg_blood", label="6-month weight change")
sigwtchgsto21all <- all_results(dir="./data_generated/maaslin_results/wtchgsto21_taxonomy/all_results.tsv",    exposure = "wtchgsto21", label="Weight change since age 21")

join1<-full_join(sigpaee_pam, sigact_paqlong, by="feature")
join2<-full_join(join1, sigbmi_dlw, by="feature")
join3<-full_join(join2, sigpfat_dlw, by="feature")
join4<-full_join(join3, sigweightchg_blood, by="feature")
join5<-full_join(join4, sigwtchgsto21, by="feature")
sigfeature<-subset(join5, select=feature)
join1<-left_join(sigfeature, sigpaee_pamall, by="feature")
join2<-left_join(sigfeature, sigact_paqlongall, by="feature")
join3<-left_join(sigfeature, sigbmi_dlwall, by="feature")
join4<-left_join(sigfeature, sigpfat_dlwall, by="feature")
join5<-left_join(sigfeature, sigweightchg_bloodall, by="feature")
join6<-left_join(sigfeature, sigwtchgsto21all, by="feature")

bind9<-rbind(join1, join2, join3, join4, join5, join6)

# read in table of species matched to phyla
phylum<-read.csv(file = './data_generated/phylum_species.csv', 
                 check.names=FALSE)
phylum<-phylum[,-c(1,4)]
colnames(phylum)[which(colnames(phylum) == 'species')] <- 'feature'
bind9<-left_join(bind9, phylum, by="feature")
bind9$feature <- substring(bind9$feature, 4)
bind9$phylum <- substring(bind9$phylum, 4)
bind9<-bind9[order(bind9$meta, bind9$phylum, bind9$feature),]
bind9$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind9$feature)
rownames(bind9) <- 1:nrow(bind9)
bug_names<-bind9[1:40, (ncol(bind9)-1):ncol(bind9)]
                                              
level_x_order <- factor(bind9$meta, level = c('Recent total PA',
                                              'Long-term total PA',
                                              'BMI',
                                              'Fat mass %',
                                              '6-month weight change',
                                              'Weight change since age 21'
                                              ))
level_y_order <- factor(bind9$name, levels = bug_names$name)
bind9$stars <- cut(bind9$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))

# create heatmap for physical activity-taxonomy associations
p1<-ggplot(bind9, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE) +
  xlab("PA, weight, and biomarkers") +
  theme(legend.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.position = "left",
        plot.title = element_text(size=50),
        axis.title.x=element_text(size=50,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=50,color="black")) +
  labs(fill = "Beta coefficient")

# create right bar for categorizing species to phyla
bind9$category<-as.factor(bind9$phylum)
p2<-ggplot(bind9, aes(level_x_order, level_y_order))+
  scale_y_discrete(position = "right")+
  geom_tile(aes(fill = category)) +
  xlab("PA, weight, and biomarkers") +
  theme(legend.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        axis.title.y= element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=50,face="bold",colour="white"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=50,colour="white")) +
  labs(fill = "Phyla")

png(file="./figure_generated/fig2c_heatmap_maaslin_taxonomy.png",width=3500,height=3500, pointsize=50)
grid.arrange(p1,p2, ncol=2, nrow=1,heights=c(8),widths=c(4, 0.7))
dev.off()


# Figure 2a: create input data for phylogenetic tree creation GraPhlan
tree<-read.delim(file='./data_generated/mlvstree.txt')
bind9$coep<-as.integer(bind9$coef*1000)

bind9coeppos <- subset(bind9, coep> 0)
bind9coepneg <- subset(bind9, coep<= 0)
bind9coepposDeciles<-quantile(bind9coeppos$coep, prob = seq(0, 1, length = 11), type = 5)
bind9coepposDeciles
#0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
#1    1    1    2    3    4    5    6    7    9   25 
bind9coepnegDeciles<-quantile(bind9coepneg$coep, prob = seq(0, 1, length = 11), type = 5)
bind9coepnegDeciles
#  0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
#-16   -5   -3   -2   -1   -1   -1    0    0    0    0 

bind9$coep <- ifelse(bind9$coep<(-6), -6, bind9$coep)
bind9$coep <- ifelse(bind9$coep>8, 8, bind9$coep)
bind9$coep <- ifelse(bind9$coep==0 &bind9$coef<0, -0.5, bind9$coep)
bind9$coep <- ifelse(bind9$coep==0 &bind9$coef>0, 0.5, bind9$coep)

bind9<-subset(bind9, select=-name)
out_rings_color <-function(meta, num){
  sig_med <- bind9[bind9$metadata==meta,]
  sig_med$species<-tolower(sig_med$feature)
  match<-full_join(tree, sig_med, by="species")
  
  match$anno2<-"ring_color"
  match$ring<-num
  match$color[match$coep==0.5]<-"#c1e0d1"
  match$color[match$coep==1]  <-"#b2d8c5"
  match$color[match$coep==2]  <-"#a3d0ba"
  match$color[match$coep==3]  <-"#93c9ae"
  match$color[match$coep==4]  <-"#84c1a3"
  match$color[match$coep==5]  <-"#75b997"
  match$color[match$coep==6]  <-"#66b28c"
  match$color[match$coep==7]  <-"#5ba07e"
  match$color[match$coep==8]  <-"#518e70"
  
  match$color[match$coep==-0.5]<-"#ffc5c5"
  match$color[match$coep==-1]  <-"#ffb2b2"
  match$color[match$coep==-2]  <-"#ff9f9f"
  match$color[match$coep==-3]  <-"#ff8c8c"
  match$color[match$coep==-4]  <-"#ff7979"
  match$color[match$coep==-5]  <-"#ff6666"
  match$color[match$coep==-6]  <-"#ff5353"
  
  match$color[is.na(match$coep)]<-"#CCCCC3"
  match$color<-toupper(match$color)
  
  color<-subset(match, select = c(name, anno2, ring, color))
  
  return(color)
}
co1<-out_rings_color(meta="paee_pam", num=7)
co2<-out_rings_color(meta = "act_paqlong" , num=6)
co3<-out_rings_color(meta = "bmi_dlw", num=5)
co4<-out_rings_color(meta = "pfat_dlw", num=4)
co5<-out_rings_color(meta = "weightchg_blood", num=3)
co6<-out_rings_color(meta = "wtchgsto21", num=2)

coall <-rbind(co1, co2, co3, co4, co5, co6)
write.table(coall, "./data_generated/out_rings_color.txt", 
            sep="\t", quote = FALSE, col.names = F, row.names = F)

out_rings_alpha <-function(path, meta, num){
  
  sig_med <- bind9[bind9$metadata==meta,]
  
  sig_med$species<-tolower(sig_med$feature)
  
  match<-full_join(tree, sig_med, by="species")
  
  match$anno1<-"ring_alpha"
  match$ring<-num
  match$alpha[is.na(match$coef)]<-0.3
  match$alpha[!is.na(match$coef)]<-0.8
  
  alpha<-subset(match, select=c(name, anno1, ring, alpha))
  return(alpha)
}
al1<-out_rings_alpha(path="./data_generated/maaslin_results/paee_pam_taxonomy/significant_results.tsv", meta="paee_pam", num=7)
al2<-out_rings_alpha(path="./data_generated/maaslin_results/act_paqlong_taxonomy/significant_results.tsv",    meta = "act_paqlong", num=6)
al3<-out_rings_alpha(path="./data_generated/maaslin_results/bmi_dlw_taxonomy/significant_results.tsv",   meta = "bmi_dlw", num=5)
al4<-out_rings_alpha(path="./data_generated/maaslin_results/pfat_dlw_taxonomy/significant_results.tsv",     meta = "pfat_dlw", num=4)
al5<-out_rings_alpha(path="./data_generated/maaslin_results/weightchg_blood_taxonomy/significant_results.tsv",  meta = "weightchg_blood", num=3)
al6<-out_rings_alpha(path="./data_generated/maaslin_results/wtchgsto21_taxonomy/significant_results.tsv",  meta = "wtchgsto21", num=2)

alall<-rbind(al1, al2, al3, al4, al5, al6)
write.table(alall, "./data_generated/out_rings_alpha.txt", 
            sep="\t", quote = FALSE, col.names = F, row.names = F)
