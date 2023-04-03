
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity and abundance of A. putredinis in relation to body weight change
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
library(lmerTest)
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
library(nlme)

meta_species<-read.csv(file="./data_generated/mbs_meta_species.csv", header = TRUE)

#s__Alistipes_putredinis
summary(meta_species$s__Alistipes_putredinis)
meta_species$s__Alistipes_putredinis2c <- with(meta_species, ifelse(s__Alistipes_putredinis> 4.5867, 1, 0))

# bmimbs #
a <- lmer(bmimbs ~ act13m + s__Alistipes_putredinis2c + act13m*s__Alistipes_putredinis2c + (1 | id) , data=meta_species)
summary(a)
#all act 
a <- lmer( bmimbs ~ act13m + s__Alistipes_putredinis2c + act13m*s__Alistipes_putredinis2c +age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id) , data=meta_species)
summary(a)
#in each stratum
a.below <- lmer ( bmimbs ~ act13m +age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(a.below)
confint(a.below)

a.above <- lmer ( bmimbs ~ act13m +age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(a.above)
confint(a.above)

meta_species$cutx[meta_species$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/mbs_act13m_bmimbs.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act13m, bmimbs, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab(expression(BMI~at~stool~collection~(kg/(m^2))))+theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_text(size=50),
        axis.title.y = element_text(size=50),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

# wtchg18mbs
a <- lmer( wtchg1813 ~ act13m + s__Alistipes_putredinis2c + act13m*s__Alistipes_putredinis2c + (1 | id) , data=meta_species)
summary(a)

a <- lmer( wtchg1813 ~ act13m + s__Alistipes_putredinis2c + act13m*s__Alistipes_putredinis2c +wt18+age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id) , data=meta_species)
summary(a)

#in each stratum
a.below <- lmer ( wtchg18mbs ~ act8913v +wt18+age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id), data=meta_species[meta_species$s__Alistipes_putredinis2c==0, ])
summary(a.below)
confint(a.below)

a.above <- lmer ( wtchg18mbs ~ act8913v +wt18+age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id), data=meta_species[meta_species$s__Alistipes_putredinis2c==1, ])
summary(a.above)
confint(a.above)

meta_species$cutx[meta_species$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species$cutx[meta_species$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/mbs_act8913v_wtchg18mbs.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species, aes(act8913v, wtchg18mbs, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("Weight change from age 21 to stool collection")+theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_text(size=50),
        axis.title.y = element_text(size=50),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()

# wtchg1317 
#in each stratum
a.below <- lmer ( wtchg1317 ~ act13m +age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id), data=meta_species2[meta_species2$s__Alistipes_putredinis2c==0, ])
summary(a.below)
confint(a.below)

a.above <- lmer ( wtchg1317 ~ act13m +age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id), data=meta_species2[meta_species2$s__Alistipes_putredinis2c==1, ])
summary(a.above)
confint(a.above)

a <- lmer( wtchg1317 ~ act13m + s__Alistipes_putredinis2c + act13m*s__Alistipes_putredinis2c +age13+calormbs+smk+probio_2mo_qu+stool_type+ant_12mo_qu+meno+(1 | id) , data=meta_species2)
summary(a)

meta_species2$cutx[meta_species2$s__Alistipes_putredinis2c==0]<-"\u2264 median"
meta_species2$cutx[meta_species2$s__Alistipes_putredinis2c==1]<-"> median"
png(file="./figure_generated/mbs_act13m_wtchg1317.png",width=1200,height=1200, pointsize=50)
ggplot(meta_species2, aes(act13m, wtchg1317, col=factor(cutx))) + geom_point(size=15, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=6) +
  ggtitle("Alistipes putredinis")+
  xlab("Physical activity level (MET-hours/week)")+
  ylab("Weight change 4 years after stool collection (kg)")+theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=50),
        legend.position="bottom",
        legend.text = element_text(size=50),
        legend.title = element_blank(),
        axis.title.x = element_text(size=50),
        axis.title.y = element_text(size=50),
        axis.text.y=element_text(size=50),
        axis.text.x=element_text(size=50))
dev.off()
