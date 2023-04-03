

*******************************************************************
* Table 2 correlations between the exposure and outcome variables *
*******************************************************************;

%inc '/udd/nhkwa/mlvs/mlvs_meta_data.sas';

ods html file = 'correlate.xls';
proc corr data = all nomiss;
var paee_pamwk mets_vigpamwk mets_modpamwk mets_ltpamwk act_paqlong vpa_paqlong mpa_paqlong lpa_paqlong bmi_dlw pfat_dlw bodyweightchg_dlw wtchgsto21 lgcrp lghba1c;
run;
ods html close;

proc means data=all n nmiss mean std min median max;
var totmets_paq meat122c fish122c wgrain122c legu122c fruit122c veg122c nut122c alco122cn calor122cn;
run;