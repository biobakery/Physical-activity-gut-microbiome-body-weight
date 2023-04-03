
*********************************************************
*                           Table 1                     *
*********************************************************;

%inc '/udd/nhkwa/mlvs/mlvs_meta_data.sas';

data all;
set all;
if bristol_cat=1 then bristol_cat1=1;else bristol_cat1=0;
if bristol_cat=2 then bristol_cat2=1;else bristol_cat2=0;
if bristol_cat=3 then bristol_cat3=1;else bristol_cat3=0;
smk=1-smk;
run;

data stool1;
set all;
if cvisit=1 then output stool1;
run;

proc rank data=all groups=4 out=all;
  var paee_pamwk;
  ranks paee_pamwkq;
run;

proc sort data=all;by paee_pamwkq;run;
proc means data=all n nmiss min median max;
var paee_pamwk;
by paee_pamwkq;
run;
proc means data=all;
var age_fecal mets_vigpamwk mets_modpamwk mets_ltpamwk act_paqlong vpa_paqlong mpa_paqlong lpa_paqlong
    bmi_dlw pfat_dlw bodyweightchg_dlw wtchgsto21 calor122cn crp_plasma hba1cp;
run;

proc freq data=all;
tables smk probio_2mo_qu ant_12mo_qu bristol_cat1 bristol_cat2 bristol_cat3;
run;

proc sort data=all NODUPKEY out=all_ex;
by SubjectID;
run; 
proc freq data=all_ex;
tables paee_pamwkq;
run;
proc means data=all_ex;
var age_fecal;
run; 

proc rank data=stool1 groups=4 out=stool1;
  var paee_pamwk;
  ranks paee_pamwkq;
run;

%table1(data=all,
        exposure=  paee_pamwkq ,
        varlist  = age_fecal mets_vigpamwk mets_modpamwk mets_ltpamwk act_paqlong vpa_paqlong mpa_paqlong lpa_paqlong
                   bmi_dlw pfat_dlw bodyweightchg_dlw wtchgsto21   
                   smk calor122cn probio_2mo_qu ant_12mo_qu bristol_cat1 bristol_cat2 bristol_cat3 
                   crp_plasma hba1cp ,
        cat      = smk probio_2mo_qu ant_12mo_qu bristol_cat1 bristol_cat2 bristol_cat3 ,
        ageadj   =  F,
        rtftitle =  Characteristics in all,
        landscape=  F,
        file     =  table1_all,
        dec      =  1,
        uselbl   =  F);
run;

%table1(data=stool1,
        exposure=  paee_pamwkq ,
        varlist  = age_fecal bodyweightchg_dlw ,
        ageadj   =  F,
        rtftitle =  Characteristics in stool1,
        landscape=  F,
        file     =  table1_stool1,
        dec      =  1,
        uselbl   =  F);
run;


