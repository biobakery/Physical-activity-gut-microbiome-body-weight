##############################################################################################################################################
# 1) Purposes: conduct MaAsLin regression to calculate the associations of the measures of physical activity and adiposity with abundances of microbial functional pathways, and create heatmap
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
library(Maaslin2)
library(ggplot2)
library(grid)
library(ggpubr)
library(haven)
library(Hmisc)
library(dplyr)
library(reshape)

# read in DNA pathway data
s1<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s1.sas7bdat")
s2<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s2.sas7bdat")
s3<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s3.sas7bdat")
s4<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s4.sas7bdat")
midnumber_aliasid<-read.csv(file="./data_generated/midnumber_aliasid.csv")
subjectID_aliasid_key<-read.csv(file="./data_generated/subjectID_aliasid_key.csv")
id_subject<-inner_join(midnumber_aliasid, subjectID_aliasid_key, by="aliasid")
colnames(id_subject)[which(names(id_subject) == "SubjectID")] <- "id"
colnames(id_subject)[which(names(id_subject) == "mid_number")] <- "SubjectID"

transform_data<-function(data1, label1, label2){
  data1<-inner_join(data1, id_subject, by="id")
  data1$sampleid<-paste0(data1$SubjectID, label1)
  data2<-as.data.frame(t(data1[,-(ncol(data1)-3):-(ncol(data1))]))
  colnames(data2)<-data1$sampleid
  rownames(data2)<-gsub(label2, "", rownames(data2))
  data2<-data2[order(rownames(data2)),]
  return(data2)
}
t_s1<-transform_data(s1, "_SF05", "_s1")
t_s2<-transform_data(s2, "_SF06", "_s2")
t_s3<-transform_data(s3, "_SF07", "_s3")
t_s4<-transform_data(s4, "_SF08", "_s4")
ec_table<-cbind(t_s1, t_s2, t_s3, t_s4)

t_relab1<-as.data.frame(t(ec_table))
t_relab1<-t_relab1[, colSums(t_relab1)!=0]
t_relab1<-t_relab1[,order(-colSums(t_relab1))]

# QC DNA pathway data: orthoganal filtering
dna.pwy.qc <- t_relab1[,c(),drop=F]

for (i in seq_len(ncol(t_relab1)))
{
  good <- T
  for (j in seq_len(ncol(dna.pwy.qc)))
  {
    if (abs(cor.test(t_relab1 [, i], dna.pwy.qc[, j], na.rm=TRUE)$estimate) > 0.9)
    {
      good <- F
    }
  }
  if (good)
  {
    dna.pwy.qc <- cbind(dna.pwy.qc, t_relab1 [, i])
    colnames(dna.pwy.qc)[ncol(dna.pwy.qc)] <-
      colnames(t_relab1)[i]
  }
}
rownames(dna.pwy.qc)[which(rownames(dna.pwy.qc)=="70324089_SF07")]<-"70324089_SF05"
rownames(dna.pwy.qc)[which(rownames(dna.pwy.qc)=="70324089_SF08")]<-"70324089_SF06"
# Only keep metagenomes with average sequencing depth greater than 1M
# the id list of high sequencing depth samples is from Abu-Ali Nat Microbiol 2018
id_913<-read.csv(file="./data_generated/id_high_depth.csv")
rownames(id_913)<-id_913$id #mid in fact
dna.pwy.qc$id<-rownames(dna.pwy.qc) #mid in fact
dna.pwy.qc<-inner_join(dna.pwy.qc, id_913, by="id")
rownames(dna.pwy.qc)<-dna.pwy.qc$id #mid in fact
dna.pwy.qc<-subset(dna.pwy.qc, select = -id)
# read in meatdata
meta_med<-read.csv(file="./data_generated/meta_species.csv")
meta_med$mvpa_paqlong <- meta_med$vpa_paqlong + meta_med$mpa_paqlong
rownames(meta_med)<-meta_med$mid
meta_med$age_fecal<-as.numeric(as.character(meta_med$age_fecal))
meta_med$calor122cn<-as.numeric(as.character(meta_med$calor122cn))
meta_med$totMETs_paq<-as.numeric(as.character(meta_med$totMETs_paq))

# maaslin regressions
maaslin_diet<-function(dir, exposure){
  Maaslin2(dna.pwy.qc, 
           meta_med, 
           dir,
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="SubjectID", 
           fixed_effects = c(exposure,"age_fecal","ahei","calor122cn","probio_2mo_qu","stool_type", "smk", "ant_12mo_qu"))
}
maaslin_diet(dir="./data_generated/maaslin_results/paee_pam_pwy_dna",  exposure = "paee_pam")
maaslin_diet(dir="./data_generated/maaslin_results/mets_vigpam_pwy_dna",  exposure = "mets_vigpam")
maaslin_diet(dir="./data_generated/maaslin_results/mets_modpam_pwy_dna",  exposure = "mets_modpam")
maaslin_diet(dir="./data_generated/maaslin_results/mets_ltpam_pwy_dna",  exposure = "mets_ltpam")
maaslin_diet(dir="./data_generated/maaslin_results/act_paqlong_pwy_dna",  exposure = "act_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/vpa_paqlong_pwy_dna",  exposure = "vpa_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/mpa_paqlong_pwy_dna",  exposure = "mpa_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/lpa_paqlong_pwy_dna",  exposure = "lpa_paqlong")
maaslin_diet(dir="./data_generated/maaslin_results/bmi_dlw_pwy_dna",  exposure = "bmi_dlw")
maaslin_diet(dir="./data_generated/maaslin_results/pfat_dlw_pwy_dna",  exposure = "pfat_dlw")
maaslin_diet(dir="./data_generated/maaslin_results/weightchg_blood_pwy_dna",  exposure = "weightchg_blood")
maaslin_diet(dir="./data_generated/maaslin_results/wtchgsto21_pwy_dna",  exposure = "wtchgsto21")
maaslin_diet(dir="./data_generated/maaslin_results/lgcrp_pwy_dna",  exposure = "lgcrp")
maaslin_diet(dir="./data_generated/maaslin_results/lghba1c_pwy_dna",  exposure = "lghba1c")

# combine resutls from maaslin regressions of different dietary variables
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}

sigpaee_pam <- sig_results(dir="./data_generated/maaslin_results/paee_pam_pwy_dna/significant_results.tsv",  exposure = "paee_pam")
sigmets_vigpam <- sig_results(dir="./data_generated/maaslin_results/mets_vigpam_pwy_dna/significant_results.tsv",  exposure = "mets_vigpam")
sigmets_modpam <- sig_results(dir="./data_generated/maaslin_results/mets_modpam_pwy_dna/significant_results.tsv",  exposure = "mets_modpam")
sigmets_ltpam <- sig_results(dir="./data_generated/maaslin_results/mets_ltpam_pwy_dna/significant_results.tsv",  exposure = "mets_ltpam")
sigact_paqlong <- sig_results(dir="./data_generated/maaslin_results/act_paqlong_pwy_dna/significant_results.tsv",  exposure = "act_paqlong")
sigvpa_paqlong <- sig_results(dir="./data_generated/maaslin_results/vpa_paqlong_pwy_dna/significant_results.tsv",  exposure = "vpa_paqlong")
sigmpa_paqlong <- sig_results(dir="./data_generated/maaslin_results/mpa_paqlong_pwy_dna/significant_results.tsv",  exposure = "mpa_paqlong")
siglpa_paqlong <- sig_results(dir="./data_generated/maaslin_results/lpa_paqlong_pwy_dna/significant_results.tsv",  exposure = "lpa_paqlong")
sigbmi_dlw <- sig_results(dir="./data_generated/maaslin_results/bmi_dlw_pwy_dna/significant_results.tsv",  exposure = "bmi_dlw")
sigpfat_dlw <- sig_results(dir="./data_generated/maaslin_results/pfat_dlw_pwy_dna/significant_results.tsv",  exposure = "pfat_dlw")
sigweightchg_blood <- sig_results(dir="./data_generated/maaslin_results/weightchg_blood_pwy_dna/significant_results.tsv",  exposure = "weightchg_blood")
sigwtchgsto21 <- sig_results(dir="./data_generated/maaslin_results/wtchgsto21_pwy_dna/significant_results.tsv",  exposure = "wtchgsto21")
siglgcrp <- sig_results(dir="./data_generated/maaslin_results/lgcrp_pwy_dna/significant_results.tsv",  exposure = "lgcrp")
siglghba1c <- sig_results(dir="./data_generated/maaslin_results/lghba1c_pwy_dna/significant_results.tsv",  exposure = "lghba1c")

all_results <- function(dir, exposure, label){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}
sigpaee_pamall <- all_results(dir="./data_generated/maaslin_results/paee_pam_pwy_dna/all_results.tsv",  exposure = "paee_pam", label="Recent total PA")
sigmets_vigpamall <- all_results(dir="./data_generated/maaslin_results/mets_vigpam_pwy_dna/all_results.tsv",  exposure = "mets_vigpam", label="Recent vigorous PA")
sigmets_modpamall <- all_results(dir="./data_generated/maaslin_results/mets_modpam_pwy_dna/all_results.tsv",  exposure = "mets_modpam", label="Recent moderate PA")
sigmets_ltpamall <- all_results(dir="./data_generated/maaslin_results/mets_ltpam_pwy_dna/all_results.tsv",  exposure = "mets_ltpam", label="Recent light PA")
sigact_paqlongall <- all_results(dir="./data_generated/maaslin_results/act_paqlong_pwy_dna/all_results.tsv",     exposure = "act_paqlong", label="Long-term total PA")
sigvpa_paqlongall <- all_results(dir="./data_generated/maaslin_results/vpa_paqlong_pwy_dna/all_results.tsv",     exposure = "vpa_paqlong", label="Long-term vigorous PA")
sigmpa_paqlongall <- all_results(dir="./data_generated/maaslin_results/mpa_paqlong_pwy_dna/all_results.tsv",     exposure = "mpa_paqlong", label="Long-term moderate PA")
siglpa_paqlongall <- all_results(dir="./data_generated/maaslin_results/lpa_paqlong_pwy_dna/all_results.tsv",     exposure = "lpa_paqlong", label="Long-term light PA")
sigbmi_dlwall <- all_results(dir="./data_generated/maaslin_results/bmi_dlw_pwy_dna/all_results.tsv",  exposure = "bmi_dlw", label="BMI at stool collection")
sigpfat_dlwall <- all_results(dir="./data_generated/maaslin_results/pfat_dlw_pwy_dna/all_results.tsv",    exposure = "pfat_dlw", label="Fat mass % at stool collection")
sigweightchg_bloodall <- all_results(dir="./data_generated/maaslin_results/weightchg_blood_pwy_dna/all_results.tsv",    exposure = "weightchg_blood", label="6-month weight change between the 2 stool collections")
sigwtchgsto21all <- all_results(dir="./data_generated/maaslin_results/wtchgsto21_pwy_dna/all_results.tsv",    exposure = "wtchgsto21", label="Weight change from age 21 to stool collection")
siglgcrpall <- all_results(dir="./data_generated/maaslin_results/lgcrp_pwy_dna/all_results.tsv",    exposure = "lgcrp", label="CRP at stool collection")
siglghba1call <- all_results(dir="./data_generated/maaslin_results/lghba1c_pwy_dna/all_results.tsv",    exposure = "lghba1c", label="HbA1c at stool collection")

join1<-full_join(sigpaee_pam, sigmets_vigpam, by="feature")
join2<-full_join(join1, sigmets_modpam, by="feature")
join3<-full_join(join2, sigmets_ltpam, by="feature")
join4<-full_join(join3, sigact_paqlong, by="feature")
join5<-full_join(join4, sigvpa_paqlong, by="feature")
join6<-full_join(join5, sigmpa_paqlong, by="feature")
join7<-full_join(join6, siglpa_paqlong, by="feature")
join8<-full_join(join7, sigbmi_dlw, by="feature")
join9<-full_join(join8, sigpfat_dlw, by="feature")
join10<-full_join(join9, sigweightchg_blood, by="feature")
join11<-full_join(join10, sigwtchgsto21, by="feature")
join12<-full_join(join11, siglgcrp, by="feature")
join13<-full_join(join12, siglghba1c, by="feature")
sigfeature<-subset(join13, select=feature)
join1<-left_join(sigfeature, sigpaee_pamall, by="feature")
join2<-left_join(sigfeature, sigmets_vigpamall, by="feature")
join3<-left_join(sigfeature, sigmets_modpamall, by="feature")
join4<-left_join(sigfeature, sigmets_ltpamall, by="feature")
join5<-left_join(sigfeature, sigact_paqlongall, by="feature")
join6<-left_join(sigfeature, sigvpa_paqlongall, by="feature")
join7<-left_join(sigfeature, sigmpa_paqlongall, by="feature")
join8<-left_join(sigfeature, siglpa_paqlongall, by="feature")
join9<-left_join(sigfeature, sigbmi_dlwall, by="feature")
join10<-left_join(sigfeature, sigpfat_dlwall, by="feature")
join11<-left_join(sigfeature, sigweightchg_bloodall, by="feature")
join12<-left_join(sigfeature, sigwtchgsto21all, by="feature")
join13<-left_join(sigfeature, siglgcrpall, by="feature")
join14<-left_join(sigfeature, siglghba1call, by="feature")

bind9<-rbind(join1, join2, join3, join4, join5, join6, join7, join8, join9, join10, join11, join12, join13, join14)

colnames(bind9)[which(colnames(bind9)=="feature")]<-"pwy_channing"
pwy_label<-read.csv(file="./data_generated/dna_pwy_channing_label.csv")
bind9_label<-left_join(bind9, pwy_label, by="pwy_channing")
write.csv(bind9_label, file="./data_generated/dna_pwy_maaslin.csv")

# create heatmap for PA-dna pathway associations
level_x_order <- factor(bind9_label$meta, level = c('Recent total PA',
                                                    'Recent vigorous PA',
                                                    'Recent moderate PA',
                                                    'Recent light PA',
                                                    'Long-term total PA',
                                                    'Long-term vigorous PA',
                                                    'Long-term moderate PA',
                                                    'Long-term light PA',
                                                    'BMI at stool collection',
                                                    'Fat mass % at stool collection',
                                                    '6-month weight change between the 2 stool collections',
                                                    'Weight change from age 21 to stool collection',
                                                    'CRP at stool collection',
                                                    'HbA1c at stool collection'))

level_y_order <- factor(bind9_label$pwy_name, levels = pwy_label$pwy_name)
bind9_label$stars <- cut(bind9_label$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))

p1<-ggplot(bind9_label, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab",midpoint = 0 ) +
  geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE) +
  xlab("PA, body weight measures, and biomarkers") +
  theme(legend.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.position = "right",
        plot.title = element_text(size=50,color="black"),
        axis.title.x=element_text(size=50,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=50,color="black")) +
  labs(fill = "Beta coefficient")
png(file="./figure_generated/sfig_heatmap_maaslin_dna_pwy.png",width=3500,height=3500, pointsize=50)
print(p1)
dev.off()

# only keep pathway for creating figure 2
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

colnames(bind9)[which(colnames(bind9)=="feature")]<-"pwy_channing"
pwy_label<-read.csv(file="./data_generated/dna_pwy_channing_label.csv")
bind9_label<-left_join(bind9, pwy_label, by="pwy_channing")

bind9display<-subset(bind9_label, pwy_name=="Adenosine ribonucleotides de novo biosynthesis"|pwy_name=="Putrescine biosynthesis Iv" |
                       pwy_name=="Inosine-5'-phosphate biosynthesis I"|pwy_name=="Urate biosynthesis/inosine 5'-phosphate degradation" | 
                       pwy_name=="Glutaryl-coa degradation"|pwy_name=="L-isoleucine biosynthesis III"|pwy_name=="Pyruvate fermentation to acetate and lactate II" )

level_x_order <- factor(bind9display$meta, level = c('Recent total PA',
                                                     'Long-term total PA',
                                                     'BMI at stool collection',
                                                     'Fat mass % at stool collection',
                                                     '6-month weight change between the 2 stool collections',
                                                     'Weight change from age 21 to stool collection'))

level_y_order <- factor(bind9display$pwy_name, levels =c("Pyruvate fermentation to acetate and lactate II",
                                                        "Glutaryl-coa degradation",
                                                        "Putrescine biosynthesis Iv",
                                                        "Adenosine ribonucleotides de novo biosynthesis",
                                                        "Inosine-5'-phosphate biosynthesis I",
                                                        "Urate biosynthesis/inosine 5'-phosphate degradation",
                                                        "L-isoleucine biosynthesis III"))
bind9display$starss <- cut(bind9display$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))

p2<-ggplot(bind9display, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=starss), color="black", size=20, show.legend = TRUE) +
  xlab("PA and body weight measures") +
  theme(legend.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.position = "right",
        plot.title = element_text(size=50,color="black",face="bold"),
        axis.title.x=element_text(size=50,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=50,color="black")) +
  labs(fill = "Beta coefficient")
png(file="./figure_generated/a.fig2b_dna_pwy.png",width=3000,height=1500, pointsize=50)
print(p2)
dev.off()


