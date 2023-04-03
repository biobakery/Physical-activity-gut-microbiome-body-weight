
##############################################################################################################################################
# 1) Purposes: plot the species composition of the EC enzymes in the 11 significant Alistipes putredinis pathways that show interactions on physical activity in relation to body weight change
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
library(RColorBrewer)

midnumber_aliasid<-read.csv(file="./data_generated/midnumber_aliasid.csv")
subjectID_aliasid_key<-read.csv(file="./data_generated/subjectID_aliasid_key.csv")
id_subject<-inner_join(midnumber_aliasid, subjectID_aliasid_key, by="aliasid")
colnames(id_subject)[which(names(id_subject) == "SubjectID")] <- "id"
colnames(id_subject)[which(names(id_subject) == "mid_number")] <- "SubjectID"
# Only keep metagenomes with average sequencing depth greater than 1M
# the id list of high sequencing depth samples is from Abu-Ali Nat Microbiol 2018
id_913<-read.csv(file="./data_generated/id_high_depth.csv")
rownames(id_913)<-id_913$id #mid in fact
putridinis_EC_list <- read.csv(file="./data_generated/putrdinis_EC_list.csv")
# read in stratified DNA enzyme data
s1<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s1.sas7bdat")
s2<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s2.sas7bdat")
s3<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s3.sas7bdat")
s4<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s4.sas7bdat")

str_ec_label<-read.csv(file="./data_generated/stratified_dna_ec_channing_label.csv")
ec_label <- read.csv(file = "./data_generated/dna_ec_channing_label.csv")
meta_med<-read.csv(file="./data_generated/meta_species.csv")
meta_med$mvpa_paqlong <- meta_med$vpa_paqlong + meta_med$mpa_paqlong
meddiet<-subset(meta_med, select = c(mid, act_paqlong))
id_913$mid<-rownames(id_913)
meddiet<-inner_join(id_913, meddiet)
meddiet<-subset(meddiet, select=-id)

transform_data<-function(data, ec, label1, label2){
  data_sub<-data[ , grepl( ec , names(data) ) ]
  data_sub$id<-data$id
  data1<-inner_join(data_sub, id_subject, by="id")
  data1$sampleid<-paste0(data1$SubjectID, label1)
  data11<-data1[,-c((ncol(data1)-3):ncol(data1))]
  rownames(data11)<-data1$sampleid
  data2<-as.data.frame(t(data11))
  rownames(data2)<-gsub(label2, "", rownames(data2))
  data2<-data2[order(rownames(data2)),]
  return(data2)
}

# creat composition figures for DNA enzymes
stack_area <- function(ec_name, name) {
  all_dna_ec1$ec_channing<-rownames(all_dna_ec1)
  all_dna_ec1<-left_join(all_dna_ec1, str_ec_label, by="ec_channing")
  rownames(all_dna_ec1)<-all_dna_ec1$ec_label
  numex<-nchar(ec_name)
  all_dna_ec1$ec_name <- as.character(all_dna_ec1$ec_name)
  all_dna_ec1$jud<-substr(all_dna_ec1$ec_label, (numex+1), (numex+1))
  all_dna_ec1<-all_dna_ec1[all_dna_ec1$jud==":",]
  all_dna_ec1<-subset(all_dna_ec1, select = -c(jud, ec_channing, ec_label,ec_name,bugs))
  all_dna_ec1<-sweep(all_dna_ec1,
                     2,
                     STATS = colSums(all_dna_ec1, na.rm = T),
                     FUN = "/")
  all_dna_ec1_noun<-all_dna_ec1[!grepl( "unclassified" , rownames(all_dna_ec1) ), ]
  mean_rel<-as.data.frame(rowSums(all_dna_ec1_noun, na.rm = T)/rowSums(all_dna_ec1_noun>=0, na.rm = T))
  #mean_rel<-as.data.frame(rowSums(all_dna_ec1_noun, na.rm = T)/rowSums(all_dna_ec1_noun==0, na.rm = T))
  mean_rel<-sweep(mean_rel,
                  2,
                  STATS = colSums(mean_rel, na.rm=T),
                  FUN = "/")
  colnames(mean_rel)<-"mean_rel"
  dim_count<-mean_rel[mean_rel>0.10]
  num<-length(dim_count)
  all_dna_ec1_noun<-cbind(all_dna_ec1_noun, mean_rel)
  all_dna_ec1_noun<-all_dna_ec1_noun[order(-all_dna_ec1_noun$mean_rel),]
  
  stratum1 <- all_dna_ec1_noun[1:num, -ncol(all_dna_ec1_noun)]
  stratum2 <-
    all_dna_ec1_noun[(num + 1):nrow(all_dna_ec1_noun), -ncol(all_dna_ec1_noun)]
  stratum3<- all_dna_ec1[grepl( "unclassified" , rownames(all_dna_ec1) ), ]
  
  rownames(stratum3)<-"s__Unclassified"
  other <- as.data.frame(t(colSums(stratum2)))
  rownames(other) <- "s__Other_species"
  stratum_cleared <- rbind(stratum1, other, stratum3)
  t_stratum <-
    as.data.frame(t(stratum_cleared))
  
  t_stratum$mid <- rownames(t_stratum)
  t_stratum1 <-inner_join(t_stratum, meddiet, by="mid")
  t_stratum1<-t_stratum1[order(t_stratum1$act_paqlong),]
  t_stratum1<-subset(t_stratum1, select=-act_paqlong)
  p.stratum <- melt(t_stratum1, id = 'mid')
  p.stratum$variable <- as.factor(p.stratum$variable)
  p.stratum$value <- as.numeric(as.character(p.stratum$value))
  
  p.stratum <- p.stratum %>% separate(variable, c(NA, "bugs"), "s__")
  p.stratum$bugs <- gsub("_", " ", p.stratum$bugs)
  p.stratum$cat<-name
  return(p.stratum)
}

t_s1_sub<-transform_data(s1, "ecd6_4_1_3", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd6_4_1_3", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd6_4_1_3", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd6_4_1_3", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
p1<-stack_area(ec_name = "6.4.1.3", name="6.4.1.3: propionyl-coa carboxylase")

t_s1_sub<-transform_data(s1, "ecd1_1_1_40", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd1_1_1_40", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd1_1_1_40", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd1_1_1_40", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
p2<-stack_area(ec_name = "1.1.1.40", name="EC 1.1.1.40: malate dehydrogenase (oxaloacetate-decarboxylating) (nadp(+))")

t_s1_sub<-transform_data(s1, "ecd2_1_1_13", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd2_1_1_13", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd2_1_1_13", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd2_1_1_13", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
p3<-stack_area(ec_name = "2.1.1.13", name="2.1.1.13: methionine synthase")

t_s1_sub<-transform_data(s1, "ecd2_3_1_31", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd2_3_1_31", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd2_3_1_31", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd2_3_1_31", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
p4<-stack_area(ec_name = "2.3.1.31", name="2.3.1.31: homoserine o-acetyltransferase")

# read in stratified RNA data for creating composition figures
s1<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s1.sas7bdat")
s2<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s2.sas7bdat")
s3<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s3.sas7bdat")
s4<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s4.sas7bdat")

str_ec_label<-read.csv(file="./data_generated/stratified_rna_ec_channing_label.csv")
rna_ec_label <- read.csv(file = "./data_generated/rna_ec_channing_label.csv")

t_s1_sub<-transform_data(s1, "ecr6_4_1_3", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecr6_4_1_3", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecr6_4_1_3", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecr6_4_1_3", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
rp1<-stack_area(ec_name = "6.4.1.3", name="6.4.1.3: propionyl-coa carboxylase")

t_s1_sub<-transform_data(s1, "ecr1_1_1_40", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecr1_1_1_40", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecr1_1_1_40", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecr1_1_1_40", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
rp2<-stack_area(ec_name = "1.1.1.40", name="EC 1.1.1.40: malate dehydrogenase (oxaloacetate-decarboxylating) (nadp(+))")

t_s1_sub<-transform_data(s1, "ecr2_1_1_13", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecr2_1_1_13", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecr2_1_1_13", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecr2_1_1_13", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
rp3<-stack_area(ec_name = "2.1.1.13", name="2.1.1.13: methionine synthase")

t_s1_sub<-transform_data(s1, "ecr2_3_1_31", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecr2_3_1_31", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecr2_3_1_31", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecr2_3_1_31", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF07")]<-"70324089_SF05"
colnames(all_dna_ec1)[which(colnames(all_dna_ec1)=="70324089_SF08")]<-"70324089_SF06"
rp4<-stack_area(ec_name = "2.3.1.31", name="2.3.1.31: homoserine o-acetyltransferase")

p111<-rbind(p1, p2,p3, p4, rp1, rp2,rp3, rp4)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs),colors=colorRampPalette(brewer.pal(8, "Accent"))(11))
color<-c('Alistipes putredinis' ='#f7022a',
         'Faecalibacterium prausnitzii' ='#7FC97F',
         'Ruminococcus bromii' ='#386CB0',
         'Other species' ='#ffffbf',
         'Unclassified'='#D3D3D3',
         'Eubacterium siraeum' ='#BEAED4',
         'Alistipes onderdonkii' ='#00FFFF',
         'Akkermansia muciniphila' ='#666666',
         'Alistipes shahii' = '#1E9DC9',
         'Odoribacter splanchnicus' ='#DA2950',
         'Methanobrevibacter smithii' ='#003EFF'
)

p111<-rbind(p1)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("6.4.1.3: propionyl-coa carboxylase"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Faecalibacterium prausnitzii"] <- "12 Faecalibacterium prausnitzii"
p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Ruminococcus bromii"] <- "13 Ruminococcus bromii"
p111$bugs[p111$bugs == "Other species"] <- "14 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "15 Unclassified"

bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(5))

color<-c('12 Faecalibacterium prausnitzii'='#7FC97F',
         '11 Alistipes putredinis'='#f7022a',
         '13 Ruminococcus bromii'='#386CB0',
         '14 Other species'='#ffffbf',
         '15 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC6.4.1.3_dna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC6.4.1.3_dna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

p111<-rbind(p2)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("1.1.1.40: malate dehydrogenase (oxaloacetate-decarboxylating) (nadp(+))"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Other species"] <- "12 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "13 Unclassified"

bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(3))

color<-c('11 Alistipes putredinis'='#f7022a',
         '12 Other species'='#ffffbf',
         '13 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC1.1.1.40_dna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC1.1.1.40_dna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

p111<-rbind(p3)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("2.1.1.13: methionine synthase"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Faecalibacterium prausnitzii"] <- "12 Faecalibacterium prausnitzii"
p111$bugs[p111$bugs == "Ruminococcus bromii"] <- "13 Ruminococcus bromii"
p111$bugs[p111$bugs == "Eubacterium siraeum"] <- "14 Eubacterium siraeum"
p111$bugs[p111$bugs == "Other species"] <- "15 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "16 Unclassified"

bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(6))

color<-c('11 Alistipes putredinis'='#f7022a',
         '12 Faecalibacterium prausnitzii'='#7FC97F',
         '13 Ruminococcus bromii'='#386CB0',
         '14 Eubacterium siraeum'='#BEAED4',
         '15 Other species'='#ffffbf',
         '16 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC2.1.1.13_dna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC2.1.1.13_dna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

p111<-rbind(p4)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("2.3.1.31: homoserine o-acetyltransferase"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Alistipes onderdonkii"] <- "12 Alistipes onderdonkii"
p111$bugs[p111$bugs == "Akkermansia muciniphila"] <- "13 Akkermansia muciniphila"
p111$bugs[p111$bugs == "Other species"] <- "14 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "15 Unclassified"

bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(5))

color<-c('11 Alistipes putredinis'='#f7022a',
         '12 Alistipes onderdonkii'='#00FFFF',
         '13 Akkermansia muciniphila'='#666666',
         '14 Other species'='#ffffbf',
         '15 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC2.3.1.31_dna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC2.3.1.31_dna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

p111<-rbind(rp1)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("6.4.1.3: propionyl-coa carboxylase"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Alistipes onderdonkii"] <- "12 Alistipes onderdonkii"
p111$bugs[p111$bugs == "Alistipes shahii"] <- "13 Alistipes shahii"
p111$bugs[p111$bugs == "Odoribacter splanchnicus"] <- "14 Odoribacter splanchnicus"
p111$bugs[p111$bugs == "Other species"] <- "15 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "16 Unclassified"

bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(6))

color<-c('11 Alistipes putredinis'='#f7022a',
         '12 Alistipes onderdonkii'='#00FFFF',
         '13 Alistipes shahii'='#1E9DC9',
         '14 Odoribacter splanchnicus'='#DA2950',
         '15 Other species'='#ffffbf',
         '16 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC6.4.1.3_rna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC6.4.1.3_rna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

p111<-rbind(rp2)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("1.1.1.40: malate dehydrogenase (oxaloacetate-decarboxylating) (nadp(+))"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Other species"] <- "12 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "13 Unclassified"

bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(3))

color<-c('11 Alistipes putredinis'='#f7022a',
         '12 Other species'='#ffffbf',
         '13 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC1.1.1.40_rna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC1.1.1.40_rna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

p111<-rbind(rp3)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("2.1.1.13: methionine synthase"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Methanobrevibacter smithii"] <- "12 Methanobrevibacter smithii"
p111$bugs[p111$bugs == "Other species"] <- "13 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "14 Unclassified"

bugs<-unique(p111$bugs)
bugs

bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(4))

color<-c('11 Alistipes putredinis'='#f7022a',
         '12 Methanobrevibacter smithii'='#003EFF',
         '13 Other species'='#ffffbf',
         '14 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC2.1.1.13_rna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC2.1.1.13_rna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

p111<-rbind(rp4)
p111$bugs<-capitalize(p111$bugs)
p111$value[p111$value=="NaN"]<-0
p111$cat_f<-factor(p111$cat, levels = c("2.3.1.31: homoserine o-acetyltransferase"))
bugs<-unique(p111$bugs)
bugs

p111$bugs[p111$bugs == "Alistipes putredinis"] <- "11 Alistipes putredinis"
p111$bugs[p111$bugs == "Alistipes onderdonkii"] <- "12 Alistipes onderdonkii"
p111$bugs[p111$bugs == "Methanobrevibacter smithii"] <- "13 Methanobrevibacter smithii"
p111$bugs[p111$bugs == "Other species"] <- "14 Other species"
p111$bugs[p111$bugs == "Unclassified"] <- "15 Unclassified"

bugs<-unique(p111$bugs)
bugs


bug_ratio_color<-data.frame(bugs=unique(p111$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(5))

color<-c('11 Alistipes putredinis'='#f7022a',
         '12 Alistipes onderdonkii'='#00FFFF',
         '13 Methanobrevibacter smithii'='#003EFF',
         '14 Other species'='#ffffbf',
         '15 Unclassified'='#d3d3d3')

png(file="./figure_generated/putre_EC2.3.1.31_rna_legend.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", face="bold", size=30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "italic"),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

png(file="./figure_generated/putre_EC2.3.1.31_rna.png",width=3000,height=1500, pointsize=50)
ggplot(p111,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color)+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 60, face = "italic"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()
