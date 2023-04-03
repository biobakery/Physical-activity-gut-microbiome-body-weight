
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity and metabolic pathways contributed by the microbial species Alistipes putredinis in relation to body weight change
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
library(lmerTest)
library(lme4)

# read in DNA pathway data
s1<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s1.sas7bdat")
s2<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s2.sas7bdat")
s3<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s3.sas7bdat")
s4<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s4.sas7bdat")

s11<-subset(s1, select = -c(id))
s11<-sweep(s11,
           2,
           STATS = rowSums(s1, na.rm = T),
           FUN = "/")
s11$id <- s1$id
s1 <- s11

s22<-subset(s2, select = -c(id))
s22<-sweep(s22,
           2,
           STATS = rowSums(s1, na.rm = T),
           FUN = "/")
s22$id <- s2$id
s2 <- s22

s33<-subset(s3, select = -c(id))
s33<-sweep(s33,
           2,
           STATS = rowSums(s1, na.rm = T),
           FUN = "/")
s33$id <- s3$id
s3 <- s33

s44<-subset(s4, select = -c(id))
s44<-sweep(s44,
           2,
           STATS = rowSums(s1, na.rm = T),
           FUN = "/")
s44$id <- s4$id
s4 <- s44

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
dna.pwy.qc$mid<-rownames(dna.pwy.qc)

# read in meatdata
meta_med<-read.csv(file="./data_generated/meta_species.csv")
meta_med$mvpa_paqlong <- meta_med$vpa_paqlong + meta_med$mpa_paqlong
rownames(meta_med)<-meta_med$mid
meta_med$age_fecal<-as.numeric(as.character(meta_med$age_fecal))
meta_med$calor122cn<-as.numeric(as.character(meta_med$calor122cn))
meta_med$totMETs_paq<-as.numeric(as.character(meta_med$totMETs_paq))

meta_med_pwy<-inner_join(dna.pwy.qc, meta_med, by="mid")
rownames(meta_med_pwy)<-meta_med_pwy$mid
meta_med_pwy$wtchg_resp2c <- with(meta_med_pwy, ifelse(meta_med_pwy$wtchg_resp<= 2, 1, 0))
table(meta_med_pwy$wtchg_resp2c)

#the 47 pways remain after qc 
#1 gpdPWY_6121 5-aminoimidazole ribonucleotide biosynthesis I
summary(meta_med_pwy$gpdPWY_6121)
meta_med_pwy$gpdPWY_61212c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6121> 3.870e-03, 1, 0))
table(meta_med_pwy$gpdPWY_61212c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61212c + act_paqlong*gpdPWY_61212c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61212c + act_paqlong*gpdPWY_61212c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#2 gpdGLYCOLYSIS glycolysis I (from glucose 6-phosphate)
summary(meta_med_pwy$gpdGLYCOLYSIS)
meta_med_pwy$gpdGLYCOLYSIS2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdGLYCOLYSIS> 0.0006742, 1, 0))
table(meta_med_pwy$gpdGLYCOLYSIS2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdGLYCOLYSIS2c + act_paqlong*gpdGLYCOLYSIS2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdGLYCOLYSIS2c + act_paqlong*gpdGLYCOLYSIS2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#3 gpdFASYN_ELONG fatty acid elongation -- saturated
summary(meta_med_pwy$gpdFASYN_ELONG)
meta_med_pwy$gpdFASYN_ELONG2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdFASYN_ELONG> 9.169e-05, 1, 0))
table(meta_med_pwy$gpdFASYN_ELONG2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdFASYN_ELONG2c + act_paqlong*gpdFASYN_ELONG2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdFASYN_ELONG2c + act_paqlong*gpdFASYN_ELONG2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#4 gpdPWY_4984 urea cycle
summary(meta_med_pwy$gpdPWY_4984)
meta_med_pwy$gpdPWY_49842c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_4984> 2.071e-04, 1, 0))
table(meta_med_pwy$gpdPWY_49842c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_49842c + act_paqlong*gpdPWY_49842c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_49842c + act_paqlong*gpdPWY_49842c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#5 gpdGLYCOGENSYNTH glycogen biosynthesis I (from ADP-D-Glucose)
summary(meta_med_pwy$gpdGLYCOGENSYNTH)
meta_med_pwy$gpdGLYCOGENSYNTH2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdGLYCOGENSYNTH> 1.261e-03, 1, 0))
table(meta_med_pwy$gpdGLYCOGENSYNTH2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdGLYCOGENSYNTH2c + act_paqlong*gpdGLYCOGENSYNTH2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdGLYCOGENSYNTH2c + act_paqlong*gpdGLYCOGENSYNTH2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#6 gpdPWY_7663 gondoate biosynthesis (anaerobic)
summary(meta_med_pwy$gpdPWY_7663)
meta_med_pwy$gpdPWY_76632c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_7663> 2.212e-03, 1, 0))
table(meta_med_pwy$gpdPWY_76632c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_76632c + act_paqlong*gpdPWY_76632c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_76632c + act_paqlong*gpdPWY_76632c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#7 gpdPWY_7221 guanosine ribonucleotides de novo biosynthesis
summary(meta_med_pwy$gpdPWY_7221)
meta_med_pwy$gpdPWY_72212c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_7221> 4.502e-03, 1, 0))
table(meta_med_pwy$gpdPWY_72212c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_72212c + act_paqlong*gpdPWY_72212c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_72212c + act_paqlong*gpdPWY_72212c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#8 gpdPYRIDNUCSYN NAD biosynthesis I (from aspartate)
summary(meta_med_pwy$gpdPYRIDNUCSYN)
meta_med_pwy$gpdPYRIDNUCSYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPYRIDNUCSYN> 1.122e-03, 1, 0))
table(meta_med_pwy$gpdPYRIDNUCSYN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPYRIDNUCSYN2c + act_paqlong*gpdPYRIDNUCSYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPYRIDNUCSYN2c + act_paqlong*gpdPYRIDNUCSYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#9 gpdPWY_5989 stearate biosynthesis II (bacteria and plants)
summary(meta_med_pwy$gpdPWY_5989)
meta_med_pwy$gpdPWY_59892c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_5989> 0.0005112, 1, 0))
table(meta_med_pwy$gpdPWY_59892c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_59892c + act_paqlong*gpdPWY_59892c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_59892c + act_paqlong*gpdPWY_59892c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#10 gpdPWY_5695 inosine 5'-phosphate degradation
summary(meta_med_pwy$gpdPWY_5695)
meta_med_pwy$gpdPWY_56952c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_5695> 2.993e-03, 1, 0))
table(meta_med_pwy$gpdPWY_56952c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_56952c + act_paqlong*gpdPWY_56952c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_56952c + act_paqlong*gpdPWY_56952c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#11 gpd1CMET2 tetrahydrofolate biosynthesis
summary(meta_med_pwy$gpd1CMET2)
meta_med_pwy$gpd1CMET22c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpd1CMET2> 2.484e-03, 1, 0))
table(meta_med_pwy$gpd1CMET22c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpd1CMET22c + act_paqlong*gpd1CMET22c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpd1CMET22c + act_paqlong*gpd1CMET22c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#12 gpdARO chorismate biosynthesis I
summary(meta_med_pwy$gpdARO)
meta_med_pwy$gpdARO2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdARO> 3.913e-03, 1, 0))
table(meta_med_pwy$gpdARO2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdARO2c + act_paqlong*gpdARO2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdARO2c + act_paqlong*gpdARO2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#13 gpdPWY_3841 folate transformations I
summary(meta_med_pwy$gpdPWY_3841)
meta_med_pwy$gpdPWY_38412c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_3841> 2.657e-03, 1, 0))
table(meta_med_pwy$gpdPWY_38412c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_38412c + act_paqlong*gpdPWY_38412c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_38412c + act_paqlong*gpdPWY_38412c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#14 gpdPEPTIDOGLYCANSYN peptidoglycan biosynthesis I (meso-diaminopimelate containing)
summary(meta_med_pwy$gpdPEPTIDOGLYCANSYN)
meta_med_pwy$gpdPEPTIDOGLYCANSYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPEPTIDOGLYCANSYN> 4.108e-03, 1, 0))
table(meta_med_pwy$gpdPEPTIDOGLYCANSYN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPEPTIDOGLYCANSYN2c + act_paqlong*gpdPEPTIDOGLYCANSYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPEPTIDOGLYCANSYN2c + act_paqlong*gpdPEPTIDOGLYCANSYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#15 gpdFAO fatty acid ??-oxidation I
summary(meta_med_pwy$gpdFAO)
meta_med_pwy$gpdFAO2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdFAO> 4.171e-05, 1, 0))
table(meta_med_pwy$gpdFAO2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdFAO2c + act_paqlong*gpdFAO2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdFAO2c + act_paqlong*gpdFAO2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#16 gpdPWY_5971 palmitate biosynthesis II (bacteria and plants)
summary(meta_med_pwy$gpdPWY_5971)
meta_med_pwy$gpdPWY_59712c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_5971> 9.482e-05, 1, 0))
table(meta_med_pwy$gpdPWY_59712c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_59712c + act_paqlong*gpdPWY_59712c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_59712c + act_paqlong*gpdPWY_59712c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#17 gpdPWY_6147 6-hydroxymethyl-dihydropterin diphosphate biosynthesis I
summary(meta_med_pwy$gpdPWY_6147)
meta_med_pwy$gpdPWY_61472c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6147> 0.0002372, 1, 0))
table(meta_med_pwy$gpdPWY_61472c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61472c + act_paqlong*gpdPWY_61472c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61472c + act_paqlong*gpdPWY_61472c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#18 gpdPWY_6700 queuosine biosynthesis
summary(meta_med_pwy$gpdPWY_6700)
meta_med_pwy$gpdPWY_67002c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6700> 3.371e-03, 1, 0))
table(meta_med_pwy$gpdPWY_67002c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_67002c + act_paqlong*gpdPWY_67002c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_67002c + act_paqlong*gpdPWY_67002c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#19 gpdDAPLYSINESYN L-lysine biosynthesis I
summary(meta_med_pwy$gpdDAPLYSINESYN)
meta_med_pwy$gpdDAPLYSINESYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdDAPLYSINESYN> 0.0005452, 1, 0))
table(meta_med_pwy$gpdDAPLYSINESYN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdDAPLYSINESYN2c + act_paqlong*gpdDAPLYSINESYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdDAPLYSINESYN2c + act_paqlong*gpdDAPLYSINESYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#20 gpdPWY_6608 guanosine nucleotides degradation III
summary(meta_med_pwy$gpdPWY_6608)
meta_med_pwy$gpdPWY_66082c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6608> 0.0003520, 1, 0))
table(meta_med_pwy$gpdPWY_66082c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_66082c + act_paqlong*gpdPWY_66082c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_66082c + act_paqlong*gpdPWY_66082c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#21 gpdPWY0_1296 purine ribonucleosides degradation
summary(meta_med_pwy$gpdPWY0_1296)
meta_med_pwy$gpdPWY0_12962c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY0_1296> 0.0028059, 1, 0))
table(meta_med_pwy$gpdPWY0_12962c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY0_12962c + act_paqlong*gpdPWY0_12962c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY0_12962c + act_paqlong*gpdPWY0_12962c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#22 gpdPANTO phosphopantothenate biosynthesis I
summary(meta_med_pwy$gpdPANTO)
meta_med_pwy$gpdPANTO2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPANTO> 3.079e-03, 1, 0))
table(meta_med_pwy$gpdPANTO2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPANTO2c + act_paqlong*gpdPANTO2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPANTO2c + act_paqlong*gpdPANTO2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#23 gpdPYRIDOXSYN pyridoxal 5'-phosphate biosynthesis I
summary(meta_med_pwy$gpdPYRIDOXSYN)
meta_med_pwy$gpdPYRIDOXSYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPYRIDOXSYN> 1.118e-03, 1, 0))
table(meta_med_pwy$gpdPYRIDOXSYN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPYRIDOXSYN2c + act_paqlong*gpdPYRIDOXSYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPYRIDOXSYN2c + act_paqlong*gpdPYRIDOXSYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#24 gpdGLUCONEO gluconeogenesis I
summary(meta_med_pwy$gpdGLUCONEO)
meta_med_pwy$gpdGLUCONEO2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdGLUCONEO> 0.0002814, 1, 0))
table(meta_med_pwy$gpdGLUCONEO2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdGLUCONEO2c + act_paqlong*gpdGLUCONEO2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdGLUCONEO2c + act_paqlong_sec*gpdGLUCONEO2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#25 gpdPWY_6282 palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)
summary(meta_med_pwy$gpdPWY_6282)
meta_med_pwy$gpdPWY_62822c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6282> 9.227e-05, 1, 0))
table(meta_med_pwy$gpdPWY_62822c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_62822c + act_paqlong*gpdPWY_62822c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_62822c + act_paqlong*gpdPWY_62822c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#26 gpdPWY_5667 CDP-diacylglycerol biosynthesis I
summary(meta_med_pwy$gpdPWY_5667)
meta_med_pwy$gpdPWY_56672c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_5667> 3.927e-03, 1, 0))
table(meta_med_pwy$gpdPWY_56672c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_56672c + act_paqlong*gpdPWY_56672c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_56672c + act_paqlong*gpdPWY_56672c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#27 gpdPWY_6123 inosine-5'-phosphate biosynthesis I
summary(meta_med_pwy$gpdPWY_6123)
meta_med_pwy$gpdPWY_61232c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6123> 1.612e-03, 1, 0))
table(meta_med_pwy$gpdPWY_61232c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61232c + act_paqlong*gpdPWY_61232c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61232c + act_paqlong*gpdPWY_61232c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#28 gpdPWY_6387 UDP-N-acetylmuramoyl-pentapeptide biosynthesis I (meso-diaminopimelate containing)
summary(meta_med_pwy$gpdPWY_6387)
meta_med_pwy$gpdPWY_63872c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6387> 4.530e-03, 1, 0))
table(meta_med_pwy$gpdPWY_63872c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_63872c + act_paqlong*gpdPWY_63872c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_63872c + act_paqlong*gpdPWY_63872c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#29 gpdUDPNAGSYN UDP-N-acetyl-D-glucosamine biosynthesis I
summary(meta_med_pwy$gpdUDPNAGSYN)
meta_med_pwy$gpdUDPNAGSYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdUDPNAGSYN> 0.0002937, 1, 0))
table(meta_med_pwy$gpdUDPNAGSYN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdUDPNAGSYN2c + act_paqlong*gpdUDPNAGSYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdUDPNAGSYN2c + act_paqlong*gpdUDPNAGSYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#30 gpdPWY_6122 5-aminoimidazole ribonucleotide biosynthesis II
summary(meta_med_pwy$gpdPWY_6122)
meta_med_pwy$gpdPWY_61222c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6122> 3.827e-03, 1, 0))
table(meta_med_pwy$gpdPWY_61222c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61222c + act_paqlong*gpdPWY_61222c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61222c + act_paqlong*gpdPWY_61222c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#31 gpdPWY_6163 chorismate biosynthesis from 3-dehydroquinate
summary(meta_med_pwy$gpdPWY_6163)
meta_med_pwy$gpdPWY_61632c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6163> 3.976e-03, 1, 0))
table(meta_med_pwy$gpdPWY_61632c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61632c + act_paqlong*gpdPWY_61632c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61632c + act_paqlong*gpdPWY_61632c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#32 gpdDTDPRHAMSYN dTDP-L-rhamnose biosynthesis I
summary(meta_med_pwy$gpdDTDPRHAMSYN)
meta_med_pwy$gpdDTDPRHAMSYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdDTDPRHAMSYN> 2.339e-03, 1, 0))
table(meta_med_pwy$gpdDTDPRHAMSYN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdDTDPRHAMSYN2c + act_paqlong*gpdDTDPRHAMSYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdDTDPRHAMSYN2c + act_paqlong*gpdDTDPRHAMSYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#33 gpdHSERMETANA L-methionine biosynthesis III
summary(meta_med_pwy$gpdHSERMETANA)
meta_med_pwy$gpdHSERMETANA2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdHSERMETANA> 6.004e-04, 1, 0))
table(meta_med_pwy$gpdHSERMETANA2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdHSERMETANA2c + act_paqlong*gpdHSERMETANA2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdHSERMETANA2c + act_paqlong*gpdHSERMETANA2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#34 gpdRIBOSYN2 flavin biosynthesis I (bacteria and plants)
summary(meta_med_pwy$gpdRIBOSYN2)
meta_med_pwy$gpdRIBOSYN22c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdRIBOSYN2> 0.0005769, 1, 0))
table(meta_med_pwy$gpdRIBOSYN22c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdRIBOSYN22c + act_paqlong*gpdRIBOSYN22c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdRIBOSYN22c + act_paqlong*gpdRIBOSYN22c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#35 gpdPWY0_1479 tRNA processing
summary(meta_med_pwy$gpdPWY0_1479)
meta_med_pwy$gpdPWY0_14792c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY0_1479> 2.337e-05, 1, 0))
table(meta_med_pwy$gpdPWY0_14792c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY0_14792c + act_paqlong*gpdPWY0_14792c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY0_14792c + act_paqlong*gpdPWY0_14792c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#36 gpdPPGPPMET ppGpp biosynthesis
summary(meta_med_pwy$gpdPPGPPMET)
meta_med_pwy$gpdPPGPPMET2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPPGPPMET> 7.050e-05, 1, 0))
table(meta_med_pwy$gpdPPGPPMET2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPPGPPMET2c + act_paqlong*gpdPPGPPMET2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPPGPPMET2c + act_paqlong*gpdPPGPPMET2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#37 gpdOANTIGEN O-antigen building blocks biosynthesis (E. coli)
summary(meta_med_pwy$gpdOANTIGEN)
meta_med_pwy$gpdOANTIGEN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdOANTIGEN> 0.0004457, 1, 0))
table(meta_med_pwy$gpdOANTIGEN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdOANTIGEN2c + act_paqlong*gpdOANTIGEN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdOANTIGEN2c + act_paqlong*gpdOANTIGEN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#38 gpdPWY_621 sucrose degradation III (sucrose invertase)
summary(meta_med_pwy$gpdPWY_621)
meta_med_pwy$gpdPWY_6212c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_621> 1.137e-03, 1, 0))
table(meta_med_pwy$gpdPWY_6212c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_6212c + act_paqlong*gpdPWY_6212c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_6212c + act_paqlong*gpdPWY_6212c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#39 gpdPWY_7199 pyrimidine deoxyribonucleosides salvage
summary(meta_med_pwy$gpdPWY_7199)
meta_med_pwy$gpdPWY_71992c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_7199> 6.697e-04, 1, 0))
table(meta_med_pwy$gpdPWY_71992c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_71992c + act_paqlong*gpdPWY_71992c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_71992c + act_paqlong*gpdPWY_71992c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#40 gpdPWY_5104 L-isoleucine biosynthesis IV
summary(meta_med_pwy$gpdPWY_5104)
meta_med_pwy$gpdPWY_51042c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_5104> 0.0005430, 1, 0))
table(meta_med_pwy$gpdPWY_51042c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_51042c + act_paqlong*gpdPWY_51042c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_51042c + act_paqlong*gpdPWY_51042c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#41 gpdNAGLIPASYN lipid IVA biosynthesis
summary(meta_med_pwy$gpdNAGLIPASYN)
meta_med_pwy$gpdNAGLIPASYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdNAGLIPASYN> 6.359e-05, 1, 0))
table(meta_med_pwy$gpdNAGLIPASYN2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdNAGLIPASYN2c + act_paqlong*gpdNAGLIPASYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdNAGLIPASYN2c + act_paqlong*gpdNAGLIPASYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#42 gpdPWY_6151 S-adenosyl-L-methionine cycle I
summary(meta_med_pwy$gpdPWY_6151)
meta_med_pwy$gpdPWY_61512c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6151> 4.128e-03, 1, 0))
table(meta_med_pwy$gpdPWY_61512c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61512c + act_paqlong*gpdPWY_61512c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_61512c + act_paqlong*gpdPWY_61512c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#43 gpdPWY_7197 pyrimidine deoxyribonucleotide phosphorylation
summary(meta_med_pwy$gpdPWY_7197)
meta_med_pwy$gpdPWY_71972c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_7197> 9.702e-04, 1, 0))
table(meta_med_pwy$gpdPWY_71972c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_71972c + act_paqlong*gpdPWY_71972c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_71972c + act_paqlong*gpdPWY_71972c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#44 gpdTRNA_CHARGING tRNA charging
summary(meta_med_pwy$gpdTRNA_CHARGING)
meta_med_pwy$gpdTRNA_CHARGING2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdTRNA_CHARGING> 0.002845, 1, 0))
table(meta_med_pwy$gpdTRNA_CHARGING2c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdTRNA_CHARGING2c + act_paqlong*gpdTRNA_CHARGING2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdTRNA_CHARGING2c + act_paqlong*gpdTRNA_CHARGING2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#45 gpdPWY_7184 pyrimidine deoxyribonucleotides de novo biosynthesis I
summary(meta_med_pwy$gpdPWY_7184)
meta_med_pwy$gpdPWY_71842c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_7184> 1.169e-03, 1, 0))
table(meta_med_pwy$gpdPWY_71842c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_71842c + act_paqlong*gpdPWY_71842c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_71842c + act_paqlong*gpdPWY_71842c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#46 gpd1CMET2 N10-formyl-tetrahydrofolate biosynthesis
summary(meta_med_pwy$gpd1CMET2)
meta_med_pwy$gpd1CMET22c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpd1CMET2> 2.484e-03, 1, 0))
table(meta_med_pwy$gpd1CMET22c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpd1CMET22c + act_paqlong*gpd1CMET22c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpd1CMET22c + act_paqlong*gpd1CMET22c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

#47 gpdPWY_7219 adenosine ribonucleotides de novo biosynthesis
summary(meta_med_pwy$gpdPWY_7219)
meta_med_pwy$gpdPWY_72192c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_7219> 6.788e-03, 1, 0))
table(meta_med_pwy$gpdPWY_72192c)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_72192c + act_paqlong*gpdPWY_72192c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
a <- lmer( wtchgsto21 ~ act_paqlong + gpdPWY_72192c + act_paqlong*gpdPWY_72192c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)


#11 modifying pwys in the 47 pways remain after qc among the 48 pwys in channing

#1 3 gpdFASYN_ELONG fatty acid elongation -- saturated
summary(meta_med_pwy$gpdFASYN_ELONG)
meta_med_pwy$gpdFASYN_ELONG2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdFASYN_ELONG> 9.169e-05, 1, 0))
table(meta_med_pwy$gpdFASYN_ELONG2c)
a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdFASYN_ELONG2c + act_paqlong_sec*gpdFASYN_ELONG2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)
confint(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdFASYN_ELONG2c + act_paqlong_sec*gpdFASYN_ELONG2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdFASYN_ELONG2c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdFASYN_ELONG2c==1, ])
summary(a)
confint(a)

#2 9 gpdPWY_5989 stearate biosynthesis II (bacteria and plants)
summary(meta_med_pwy$gpdPWY_5989)
meta_med_pwy$gpdPWY_59892c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_5989> 0.0005112, 1, 0))
table(meta_med_pwy$gpdPWY_59892c)
a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_59892c + act_paqlong_sec*gpdPWY_59892c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_59892c + act_paqlong_sec*gpdPWY_59892c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_59892c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_59892c==1, ])
summary(a)
confint(a)

#3 13 gpdPWY_3841 folate transformations I
summary(meta_med_pwy$gpdPWY_3841)
meta_med_pwy$gpdPWY_38412c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_3841> 2.657e-03, 1, 0))
table(meta_med_pwy$gpdPWY_38412c)
a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_38412c + act_paqlong_sec*gpdPWY_38412c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_38412c + act_paqlong_sec*gpdPWY_38412c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_38412c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_38412c==1, ])
summary(a)
confint(a)

#4 15 gpdFAO fatty acid ??-oxidation I
summary(meta_med_pwy$gpdFAO)
meta_med_pwy$gpdFAO2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdFAO> 4.171e-05, 1, 0))
table(meta_med_pwy$gpdFAO2c)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdFAO2c + act_paqlong_sec*gpdFAO2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdFAO2c + act_paqlong_sec*gpdFAO2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdFAO2c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdFAO2c==1, ])
summary(a)
confint(a)

#5 22 gpdPANTO phosphopantothenate biosynthesis I
summary(meta_med_pwy$gpdPANTO)
meta_med_pwy$gpdPANTO2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPANTO> 3.079e-03, 1, 0))
table(meta_med_pwy$gpdPANTO2c)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPANTO2c + act_paqlong_sec*gpdPANTO2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPANTO2c + act_paqlong_sec*gpdPANTO2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPANTO2c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPANTO2c==1, ])
summary(a)
confint(a)

#6 24 gpdGLUCONEO gluconeogenesis I
summary(meta_med_pwy$gpdGLUCONEO)
meta_med_pwy$gpdGLUCONEO2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdGLUCONEO> 0.0002814, 1, 0))
table(meta_med_pwy$gpdGLUCONEO2c)

a <- lmer( wtchgsto21 ~ act_paqlong + gpdGLUCONEO2c + act_paqlong*gpdGLUCONEO2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdGLUCONEO2c + act_paqlong_sec*gpdGLUCONEO2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdGLUCONEO2c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdGLUCONEO2c==1, ])
summary(a)
confint(a)

#7 25 gpdPWY_6282 palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)
summary(meta_med_pwy$gpdPWY_6282)
meta_med_pwy$gpdPWY_62822c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_6282> 9.227e-05, 1, 0))
table(meta_med_pwy$gpdPWY_62822c)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_62822c + act_paqlong_sec*gpdPWY_62822c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_62822c + act_paqlong_sec*gpdPWY_62822c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_62822c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_62822c==1, ])
summary(a)
confint(a)

#8 29 gpdUDPNAGSYN UDP-N-acetyl-D-glucosamine biosynthesis I
summary(meta_med_pwy$gpdUDPNAGSYN)
meta_med_pwy$gpdUDPNAGSYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdUDPNAGSYN> 0.0002937, 1, 0))
table(meta_med_pwy$gpdUDPNAGSYN2c)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdUDPNAGSYN2c + act_paqlong_sec*gpdUDPNAGSYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdUDPNAGSYN2c + act_paqlong_sec*gpdUDPNAGSYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdUDPNAGSYN2c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdUDPNAGSYN2c==1, ])
summary(a)
confint(a)

#9 32 gpdDTDPRHAMSYN dTDP-L-rhamnose biosynthesis I
summary(meta_med_pwy$gpdDTDPRHAMSYN)
meta_med_pwy$gpdDTDPRHAMSYN2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdDTDPRHAMSYN> 2.339e-03, 1, 0))
table(meta_med_pwy$gpdDTDPRHAMSYN2c)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdDTDPRHAMSYN2c + act_paqlong_sec*gpdDTDPRHAMSYN2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdDTDPRHAMSYN2c + act_paqlong_sec*gpdDTDPRHAMSYN2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdDTDPRHAMSYN2c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdDTDPRHAMSYN2c==1, ])
summary(a)
confint(a)

#10 44 gpdTRNA_CHARGING tRNA charging
summary(meta_med_pwy$gpdTRNA_CHARGING)
meta_med_pwy$gpdTRNA_CHARGING2c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdTRNA_CHARGING> 0.002845, 1, 0))
table(meta_med_pwy$gpdTRNA_CHARGING2c)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdTRNA_CHARGING2c + act_paqlong_sec*gpdTRNA_CHARGING2c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdTRNA_CHARGING2c + act_paqlong_sec*gpdTRNA_CHARGING2c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdTRNA_CHARGING2c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdTRNA_CHARGING2c==1, ])
summary(a)
confint(a)

#11 45 gpdPWY_7184 pyrimidine deoxyribonucleotides de novo biosynthesis I
summary(meta_med_pwy$gpdPWY_7184)
meta_med_pwy$gpdPWY_71842c <- with(meta_med_pwy, ifelse(meta_med_pwy$gpdPWY_7184> 1.169e-03, 1, 0))
table(meta_med_pwy$gpdPWY_71842c)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_71842c + act_paqlong_sec*gpdPWY_71842c + (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer( wtchgsto21 ~ act_paqlong_sec + gpdPWY_71842c + act_paqlong_sec*gpdPWY_71842c +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+ (1 | SubjectID) , data=meta_med_pwy)
summary(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_71842c==0, ])
summary(a)
confint(a)

a <- lmer ( wtchgsto21 ~ act_paqlong_sec +age_fecal+calor122cnr+smk+probio_2mo_qu+stool_type+ant_12mo_qu+ahei+(1|SubjectID), data=meta_med_pwy[meta_med_pwy$gpdPWY_71842c==1, ])
summary(a)
confint(a)
