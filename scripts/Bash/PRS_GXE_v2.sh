


mkdir -p $rawLoc
mkdir -p $scratchLoc
mkdir -p $resultsLoc
mkdir -p $mtagscript

##############################################################################
# (II) Extract Conventional predictors:
#######################################
# create a backup for reproducibility, in case it gets updated
# cp /rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/20230522/Stata/analysis.dta /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/analysis.dta 

R
library(readstata13)
library('data.table')
ceu <- read.dta13("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/analysis.dta")
setDT(ceu)
vars = attributes(ceu)$var.labels
names(vars) = names(ceu)
myData=ceu[,.(eid=idno, hypdbin)]
#write.table(vars, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/CEU_covars.txt", row.names = T, col.names = T, quote = FALSE)

write.table(vars, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/CEU_covars.txt", row.names = T, col.names = T, quote = T, sep="\t")


# write eids into a separate file
myData=ceu[,.(eid=idno)]
write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_ids", row.names = F, col.names = F, quote = F, sep="\t")



numericBinary = NULL
numericBinary_key = NULL # this has the signature of "variable name"

factors = NULL
factors_key = NULL # this has the signature of "variable name" "'factor1','factor2',...,'factorn'"
i=15

for (i in 12:227 ) {
print(paste0(i,": ", colnames(ceu)[i]) )
}


#options(error=traceback)
#options(warn=2)
for (i in 12:227 ) { # 227# 12 is the age, and 227 is the wage, the rest are non-interesting variables
print(i)
values = ceu[[ colnames(ceu)[i] ]]
treatAsNumeric=F
if( i >= 106 && i <= 108) {next} # we skip the dates, as these would get converted incorrectly, and we don't need them, as we have these converted to ages as the next few vars
if( i >= 192 && i <= 193) {next}  # also exclude subjecteplist  and subjecteplistf
if( i == 20 || i == 39 || i == 105) {next}  # also exclude married, hxochd, survtype   , as it is all missing

# if it only has missing for more than half then skip
missingVals =  is.na(values) 
if(sum(missingVals) >= length(missingVals) ) {print(paste0( colnames(ceu)[i], " has too many missing values, we skip it") ); next; }



if ( is.factor( ceu[[ colnames(ceu)[i] ]] ) ) {
print(paste0(i, ": ", colnames(ceu)[i], " is a factor")) # vars[i]
factorLevels = levels( ceu[[ colnames(ceu)[i] ]] )

# if its a factor, we find missing values, and set them to the mode factor
missingVals =  is.na(values) 
if( sum(missingVals) > 0 ) {
table_values = table(values)
mostCommon = names(table_values)[which( table_values== max(table_values))]
values[missingVals] = mostCommon

}

# if it has only 2 factor levels, then we can still code it as numeric
if(length(factorLevels) <= 2) { treatAsNumeric =T}  
} else { # if it was truly numeric to begin with
treatAsNumeric = T

# find mean and SD for filtering
missingVals =  is.na(values) 
meanValues = mean(values, na.rm = T)
sdValues = sd(values, na.rm = T)
# missing values are mean imputed
if( sum(missingVals) > 0 ) {
values[missingVals] = meanValues
}

# extreme values (3DS > mean) are set to the mean (ie missing)
extremeVals = which( abs(values) > (meanValues + 3 * sdValues) )
if( sum(extremeVals) > 0 ) {
values[extremeVals] = meanValues
}

}


# convert everything to numeric representation (for numerics leave as is, but for factors, we recode them as 0,1,2 etc)
# prepare key files:
if(treatAsNumeric) {

numericBinary = cbind(numericBinary, values)
numericBinary_key = rbind(numericBinary_key,  colnames(ceu)[i])

} else { # factors get saved into a separate DF
# want to drop unused factor levels (eg for nation, we have 213 levels, but only use 11 )
values <- factor(values)
factors = cbind(factors,as.numeric(values))
factorLevels = levels( values )


newEntry = cbind(colnames(ceu)[i], paste(shQuote(factorLevels), collapse=","))
factors_key = rbind(factors_key,newEntry )
}


#print(paste0(print(colnames(ceu)[i]), " treatAsNumeric: ", treatAsNumeric)) # vars[i]

}

# write to disk
write.table(numericBinary, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary", sep="\t", row.names = F, col.names = F, quote = FALSE)
write.table(numericBinary_key, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary_key", sep="\t", row.names = F, col.names = F, quote = FALSE)
 
write.table(factors, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors", sep="\t", row.names = F, col.names = F, quote = FALSE)
write.table(factors_key, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors_key", sep="\t", row.names = F, col.names = F, quote = FALSE)
 
quit()


wc -l $COVS_numericBinary_key # 91 numeric or binary covariates
wc -l $COVS_factors_key # 115 multi level factor covariates


###################################




################################
# 2. run PGSCatalogREST_v3.R offline to get PGS that have 50K training samples and exclude the UKB
# Manually process this to match it against ICD/self reported codes ( use the code_labels.txt/ICD10_codes.tsv/ICD9_codes.tsv/OPCS4_codes.tsv/OPCS3_codes.tsv to match the phenotypes)
# 3 types of phenos will be found: PGSlist_50K_subset_annotated.csv
# 1) with codes that could be used with the 'endpoint' script
# 2) those where a curated phenp exists aready among the covariates
# 3) not available
# use endpointExporter.R of the  phenos in category (1) to generate endpoint definition files and upload to:
# /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/endpoints/

################################
# 3. Get target phenotypes from UKBB

# self report codes:
# /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/data/curated/ukbiobank/self_report
# ICD/OPS codes: /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/data/curated/ukbiobank/hospital_records/ICD10_codes.tsv


# extraction scripts only work from here
# this is now moved to /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/deprecated/endpoints/output


# we dont include (2) as we extract that from the covariates file or manually from instructions from Scott
cd $endpointLoc
pheno_name="EFO_0003877"
arraylength=${#phearray_noT2D[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_noT2D[$j-1]}

definition=$endpointLoc$'endpoints/'$pheno_name'/endpoint_definition.txt'
# echo $definition

outfile=$rawLoc$pheno_name$'_all'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
echo "DOES NOT EXIST: "$outfile
export_UKBB_pheno_nodupes3 $pheno_name $definition $rawLoc
fi

done # end of pheno extract loop



################
# also manually  extract do those that have the phenos recorded in the CEU covars file, and export it in the same format as the rest of the phenos

R
IDs = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_ids", sep="\t")
numericBinary = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary", sep="\t")
numericBinary_key = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary_key", sep="\t")


#Alcohol consumption					EFO_0007878	behavioural	alcallbin
#Smoking Status							EFO_0006527	smoking	smallbin
#Cigarettes smoked per day				EFO_0006525	smoking	smallamt
#Height									EFO_0004339	biometric	ht
#Fasting glucose						EFO_0004465	Cardio-metabolic	glucose1
#Glycated haemoglobin levels (HbA1c)	EFO_0004541	Cardio-metabolic	hba1c

findExtremes = function(phe, limit =3) {
sds = sd(phe, na.rm = T)
means  =mean(phe, na.rm = T) 
extreIndices = which(phe > means+ sds * limit | phe <   means - sds * limit )
phe[extreIndices]
return(extreIndices)
}


GenerateResults = function(phe, pheName, IDs, limit =3) {
extreIndices = findExtremes(phe,limit) 
if(length(extreIndices) > 0) {

# z-score WITHOUT the outliers
sdWithoutExtremes = sd( phe[-extreIndices] , na.rm = T)
meanWithoutExtremes = mean( phe[-extreIndices] , na.rm = T)
phe = (phe-meanWithoutExtremes) / sdWithoutExtremes
phe[extreIndices] = -9 # set the outliers to be missing via special PLINK code '-9'
print(paste0("excluding ", length(extreIndices), " extreme phenos"))
} else {
sdWithoutExtremes = sd( phe , na.rm = T)
meanWithoutExtremes = mean( phe, na.rm = T )
phe = (phe-meanWithoutExtremes) / sdWithoutExtremes
}

whereNAsAre = is.na(phe)
phe[whereNAsAre] = -9 # also set missing data to the special PLINK code '-9'

# write out SD/Mean in case we need them later 
write.table(cbind(meanWithoutExtremes, sdWithoutExtremes),paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/",pheName, "_meanSD"), sep="\t", row.names = F, col.names = F, quote = FALSE)
phenoData = cbind(IDs,IDs,phe)
write.table(phenoData, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/",pheName), sep="\t", row.names = F, col.names = F, quote = FALSE)
}


#Height									EFO_0004339	biometric	ht
phe_index = which(numericBinary_key[,1] == "ht" )
phe = numericBinary[,phe_index]
GenerateResults(phe, "EFO_0004339_all", IDs, 4) # for height we only exclude those that are 4 sds 

#Cigarettes smoked per day				EFO_0006525	smoking	smallamt
phe_index = which(numericBinary_key[,1] == "smallamt" )
phe = numericBinary[,phe_index]
GenerateResults(phe, "EFO_0006525_all", IDs, 99) # for cigarettes smoked we don't exclude, as people who  dont smoke will be 0



#Alcohol consumption					EFO_0007878	behavioural	alcallbin
phe_index = which(numericBinary_key[,1] == "alcallbin" )
phe = numericBinary[,phe_index]
whereNAsAre = is.na(phe)
phe[whereNAsAre] = -9 # also set missing data to the special PLINK code '-9'
phenoData = cbind(IDs,IDs,phe)
write.table(phenoData, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/","EFO_0007878_all"), sep="\t", row.names = F, col.names = F, quote = FALSE)




#Smoking Status							EFO_0006527	smoking	smallbin
phe_index = which(numericBinary_key[,1] == "smallbin" )
phe = numericBinary[,phe_index]
whereNAsAre = is.na(phe)
phe[whereNAsAre] = -9 # also set missing data to the special PLINK code '-9'
phenoData = cbind(IDs,IDs,phe)
write.table(phenoData, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/","EFO_0006527_all"), sep="\t", row.names = F, col.names = F, quote = FALSE)



########################

#Fasting glucose						EFO_0004465	Cardio-metabolic	glucose1
#phe_index = which(numericBinary_key[,1] == "glucose1" )
#phe = numericBinary[,phe_index]
#GenerateResults(phe, "EFO_0004465", IDs, 4) # 

#Glycated haemoglobin levels (HbA1c)	EFO_0004541	Cardio-metabolic	hba1c
#phe_index = which(numericBinary_key[,1] == "hba1c" )
#phe = numericBinary[,phe_index]
#GenerateResults(phe, "EFO_0004541", IDs, 4) # 
########################


# instead of doing the above, we use Scott's data
cigsScott =  read.table("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/biomarkers/output/biomarkers.txt", sep="\t", header=T)
head(cigsScott)


AveragePhenos = function(cigsScott, targetPhe) {
	targetPhe_idx = which(colnames(cigsScott) == targetPhe )
	eids = c()
	phenos = c() 
	numObservarionts = c() # the per eid observation
	i=3
	for (i in 1:nrow(cigsScott) ) {
		if( i %% 5000 == 1 ) {print(paste0(i,"/",nrow(cigsScott)))}

		idx = match(cigsScott$eid[i], eids ) # find out if we have had this eid before
		currentValue= cigsScott[i,targetPhe_idx]
		incrementAmount = 1
		if (is.na(currentValue) ) incrementAmount = 0;

		if(is.na(idx) == F)  { # if we already have this eid
			if(incrementAmount != 0) { # only do anything if new value is not NA
			# handle edge-case, when the previous/first value was an NA too
			if( is.na(phenos[idx]) ) { 
			phenos[idx] = currentValue # we overwrite with the current value (as otherwise if we add anything to NA, we get NA...
			} else { # if the value before was NOT NA
			phenos[idx] = phenos[idx] + currentValue 
			} # we can add safely
			
			numObservarionts[idx] = numObservarionts[idx] +incrementAmount # increment the number of times we have seen it
			}
		} else { # first time seeing this eid
			eids = c(eids,cigsScott$eid[i] )
			phenos = c(phenos, currentValue) # this could be NA
			numObservarionts = c(numObservarionts, incrementAmount) #  we start with 1, if its not NA
		}
	}

	# now take average
	phenos = phenos/numObservarionts

	return (cbind(eids, phenos) )
}



targetPhe = "fasting_glucose"
fasting_glucose = AveragePhenos(cigsScott, targetPhe)
phe = fasting_glucose[,2]
GenerateResults(phe, "EFO_0004465_all", fasting_glucose[,1], 4)


# do we want to merge this against the CEU covars? (so that they are aligned? but does that matter, if there are a different number of indis?
hba1c = AveragePhenos(cigsScott, "hba1c")
phe = hba1c[,2]
head(phe) # 47.2 40.1 31.6 39.6 33.2 36.9  ... so use this is instead of the CEU, asthis is the 'natural' scale
GenerateResults(phe, "EFO_0004541_all", hba1c[,1], 4)

# Use Scott's definition of T2D, as he got more cases via checking incidents and also used the more sophisticated  Eastwood 2016 algorithms (https://pubmed.ncbi.nlm.nih.gov/27631769/)
t2dScott =  read.table("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/Eastwood_diabetes/output/prevalent_diabetes.txt", sep="\t", header=T)
head(t2dScott)
table(t2dScott$adjudicated_diabetes)
           # Diabetes unlikely Possible gestational diabetes
                       # 546035                           878
     # Possible type 1 diabetes      Possible type 2 diabetes
                          # 426                          3550
     # Probable type 1 diabetes      Probable type 2 diabetes
                         # 1679                         23733

t2dScott_incident =  read.table("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/Eastwood_diabetes/output/incident_diabetes.txt", sep="\t", header=T)
table(t2dScott_incident$adjudicated_diabetes)
# Incident diabetes of uncertain type             No evidence of diabetes
                               # 1878                              522038
  # Possible incident type 1 diabetes   Possible incident type 2 diabetes
                                # 260                                6718
                 # Prevalent diabetes   Probable incident type 1 diabetes
                              # 29048                                 295
  # Probable incident type 2 diabetes
                              # 16064
							  
							  
# cases: Possible type 2 diabetes + Probable type 2 diabetes
# cases: Possible incident type 2 diabetes + Probable incident type 2 diabetes
# controls: everyone else

prevalent_idx = which( t2dScott$adjudicated_diabetes == "Possible type 2 diabetes" | t2dScott$adjudicated_diabetes == "Probable type 2 diabetes")
length(prevalent_idx) #  27283
prevalents = t2dScott$eid[prevalent_idx]

incident_idx = which( t2dScott_incident$adjudicated_diabetes == "Possible incident type 2 diabetes" | t2dScott_incident$adjudicated_diabetes == "Probable incident type 2 diabetes")
length(incident_idx) # 22782
incidents = t2dScott_incident$eid[incident_idx]

all_t2dcases = c(prevalents, incidents)
length(all_t2dcases)  # 50065

#elimnate the duplicates
all_t2dcases_unique = unique(all_t2dcases)
length(all_t2dcases_unique) # 45621

# define controls, all the indis who are NOT on the above list
all_t2_controls = setdiff(t2dScott$eid, all_t2dcases)
length(t2dScott$eid)  # 576301
length(all_t2_controls) # 456839
               
#elimnate the duplicates
all_t2_controls_unique = unique(all_t2_controls)
length(all_t2_controls_unique) #  456839


t2df_cases = cbind(all_t2dcases_unique, all_t2dcases_unique, 2)
t2df_controls = cbind(all_t2_controls_unique, all_t2_controls_unique, 1)

t2df = rbind(t2df_cases,t2df_controls)
head(t2df)
nrow(t2df) # 502460

write.table(t2df, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/","MONDO_0005148_all"), sep="\t", row.names = F, col.names = F, quote = FALSE)


quit()






################################




# Exclusion lists: exclude covariates that could potentially overlap with any of the phenos I'll be training on
# factors:
# /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/CEU_covars.txt


# exclude these from all!
# repno Repeat number
# surno Survey number
# survtype Survey report type
# dos Date of survey
# doo1 Date of event/censoring
# doof Date of event/censoring (all-cause mortality)
# outage1 Age at event/censoring (yrs)
# outagef Age at event/censoring (yrs) (all-cause mortality)
# duration1 Time to event/censoring (yrs)
# durationf Time to event/censoring (yrs) (all-cause mortality)
# anyfirstevent Any-first-event indicator
# anyfatalevent Any-fatal-event indicator




# cancer-related
# epf_brecan_f Breast cancer (all-cause mortality)
# epf_digcan_f Digestive related cancer (all-cause mortality)
# epf_gencan_f Genitourinary related cancer (all-cause mortality)
# ep1_tumo_f All tumour (fatal)
# ep1_digcan_f Digestive related cancer (fatal)
# ep1_luncan_f Lung cancer (fatal)
# ep1_gencan_f Genitourinary related cancer (fatal)
# ep1_brecan_f Breast cancer (fatal)
# epf_luncan_f Lung cancer (all-cause mortality)

# height-related
# whtr Waist/height ratio
# ht Height (cm)
# bmi BMI (kg/m2)


# smoking-related
# smcigstat Smoking status: cigarettes
# smpigastat Smoking status: pipes & cigars
# smunkstat Smoking status: unknown type
# smallstat Smoking status: combined
# smallbin Smoking status
# smcigpkyr Smoking amount: cigarettes (pack years)
# smpigapkyr Smoking amount: pipes & cigars (pack years)
# smunkpkyr Smoking amount: unknown type (pack years)
# smallpkyr Smoking amount: combined (pack years)
# smallamt Smoking amount: combined (cigarettes/day)


# acholol-related
# alcallstat Alcohol status: combined
# alcallbin Alcohol status
# alcallamt Alcohol amount: combined
# alcallfreq Alcohol frequency: combined (days/week)


#  cvd related, ie Venous thromboembolism
# ep1_cv All cardiovascular
# ep1_cv_f All cardiovascular (fatal)
# ep1_cv_nf All cardiovascular (non-fatal)
# hxohd History of other HD
# lipdstat Drug status: lipid-lowering unspecified
# lipdbin Lipid-lowering unspecified drug status



# stroke-related (this includes cvd), TIA= transient ischaemic attack (TIA) or ischaemic stroke 
# hxstroke History of stroke
# hxstisc History of ischaemic stroke
# hxsthae History of haemorrhagic stroke
# hxtia History of TIA
# ep1_stri Ischaemic stroke
# ep1_stri_nf Ischaemic stroke (non-fatal)
# ep1_stri_f Ischaemic stroke (fatal)
# ep1_strh Haemorrhagic stroke
# ep1_strh_f Haemorrhagic stroke (fatal)
# ep1_strh_nf Haemorrhagic stroke (non-fatal)
# ep1_stru Unclassified stroke
# ep1_stru_f Unclassified stroke (fatal)
# ep1_stru_nf Unclassified stroke (non-fatal)
# epf_stri_f Ischaemic stroke (all-cause mortality)
# epf_strh_f Haemorrhagic stroke (all-cause mortality)
####### hxstroke_p Family history of stroke - parents # NO, keep it, same as T2D


# diabetes related
# hxdiab History of diabetes
# hxdiabbin History of diabetes
# diadstat Drug status: anti-diabetics # are we sure??? ie if someone is taking anti-diabetic meds, then they probably have diabetes
##### hxdiab_p Family history of diabetes - parents # if their parent had it, then this is NOT exclusionary




# hypertension-related
# hxhyp History of hypertension
# hypdstat Drug status: anti-hypertensives
# hypdbin Anti-hypertensives drug status


# Fasting glucose related
# glucose1 Glucose (mmol/l)

# Glycated haemoglobin levels (HbA1c)
# hba1c HbA1c (%)

###################################
R

numericBinary = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary", sep="\t")
numericBinary_key = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary_key", sep="\t")
factors = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors", sep="\t")
factors_key = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors_key", sep="\t")
# need them as list to be able to match
numericBinary_key_list = unlist( numericBinary_key )
factors_key_list = unlist(factors_key)

#### all_exclusions, "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx", "xxx"
all_exclusions = c("repno","surno","survtype","dos","doo1","doof","outage1","outagef","duration1","durationf","anyfirstevent","anyfatalevent" )
cancer_exclusions  = c(all_exclusions, "epf_brecan_f", "epf_digcan_f", "epf_gencan_f", "ep1_tumo_f", "ep1_digcan_f", "ep1_luncan_f", "ep1_gencan_f", "ep1_brecan_f", "epf_luncan_f")
height_exclusions = c(all_exclusions, "whtr", "ht", "bmi")
smoking_exclusions = c(all_exclusions, "smcigstat", "smpigastat", "smunkstat", "smallstat", "smallbin", "smcigpkyr", "smpigapkyr", "smunkpkyr", "smallpkyr", "smallamt")
alcohol_exclusions = c(all_exclusions, "alcallstat", "alcallbin", "alcallamt", "alcallfreq")
cvd_exclusions = c(all_exclusions, "ep1_cv", "ep1_cv_f", "ep1_cv_nf", "hxohd", "lipdstat", "lipdbin")
stroke_exclusions = c(cvd_exclusions, "hxstroke", "hxstisc", "hxsthae", "hxtia", "ep1_stri", "ep1_stri_nf", "ep1_stri_f", "ep1_strh", "ep1_strh_f", "ep1_strh_nf", "ep1_stru", "ep1_stru_f", "ep1_stru_nf", "epf_stri_f", "epf_strh_f")
t2d_exclusions = c(all_exclusions, "hxdiab", "hxdiabbin", "diadstat")
hypertension_exclusions = c(cvd_exclusions, "hxhyp", "hypdstat", "hypdbin")
glucose_exclusions = c(all_exclusions, "glucose1")
hba1c_exclusions = c(all_exclusions, "hba1c")

allExclusions = list (cancer_exclusions, hypertension_exclusions, stroke_exclusions, all_exclusions, cvd_exclusions, t2d_exclusions, alcohol_exclusions, smoking_exclusions, height_exclusions, glucose_exclusions, hba1c_exclusions )
covariateExclusionList=c( 'cancer', 'hypertension', 'stroke', 'all', 'cvd', 'diabetes', 'alcohol', 'smoking', 'height', 'glucose', 'hba1c' )



excludeItems = function(current_exclusions, exclType) {
# find the indices, sort them into factor/numeric categories, depending on which list they are on
factor_exclusions_indices = c()
numeric_exclusions_indices = c()

item = "urea"
i=1
for (i in 1:length(current_exclusions)) {
item = current_exclusions[i]
print(item)

indexNumeric = match(item, numericBinary_key_list )
indexFactor = match(item,factors_key_list)

if(is.na(indexNumeric) == F) { numeric_exclusions_indices = c(numeric_exclusions_indices, indexNumeric ) }
if(is.na(indexFactor) == F) { factor_exclusions_indices = c(factor_exclusions_indices, indexFactor ) }
}

# check if both lists are non-0 length
# then exclude the exclusion lists
if (length(numeric_exclusions_indices) > 0 ) {
numericBinary_exclusions = numericBinary[,-numeric_exclusions_indices]
numericBinary_key_exclusions = numericBinary_key[-numeric_exclusions_indices,]
} else {
numericBinary_exclusions = numericBinary
numericBinary_key_exclusions = numericBinary_key
}

if (length(factor_exclusions_indices) > 0 ) {
factors_exclusions = factors[,-factor_exclusions_indices]
factors_key_exclusions = factors_key[-factor_exclusions_indices,]
} else {
factors_exclusions = factors
factors_key_exclusions = factors_key
}

 
# write to disk
write.table(numericBinary_exclusions, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary_", exclType), sep="\t", row.names = F, col.names = F, quote = FALSE)
write.table(numericBinary_key_exclusions, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary_key_", exclType ), sep="\t", row.names = F, col.names = F, quote = FALSE)
 
write.table(factors_exclusions, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors_", exclType), sep="\t", row.names = F, col.names = F, quote = FALSE)
write.table(factors_key_exclusions, paste0("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors_key_", exclType), sep="\t", row.names = F, col.names = F, quote = FALSE)
 
}


j=1

for (j in 1:length(covariateExclusionList) ) {
current_exclusions = allExclusions[[j]]
exclType= covariateExclusionList[j]

excludeItems(current_exclusions,exclType)
}


####
# generate text for a human readable names for these

vars = read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/CEU_covars.txt", skip = 1, sep = "\t", quote = "\"")

factors_key_filtered = as.character(unlist(factors_key$V1)) 
factors_key_filtered = factors_key_filtered[!(factors_key_filtered %in% all_exclusions)]
length(factors_key_filtered)

numericBinary_key_filtered = as.character(unlist(numericBinary_key$V1))
numericBinary_key_filtered = numericBinary_key_filtered[!(numericBinary_key_filtered %in% all_exclusions)]
length(numericBinary_key_filtered) # 83


genText = function(currentKey) {
numericsText = ""
buffText=""
for (i in 1:length(currentKey) ) {
idx = which( vars$V1 == currentKey[i])
if(i != 1) {buffText=', '}
numericsText = paste0(numericsText, buffText ,vars$V2[idx] )
}
numericsText
return(numericsText)
}
 
genText(numericBinary_key_filtered)
genText(factors_key_filtered)

cancer_exclusions = cancer_exclusions[!(cancer_exclusions %in% all_exclusions)]
height_exclusions = height_exclusions[!(height_exclusions %in% all_exclusions)]
smoking_exclusions = smoking_exclusions[!(smoking_exclusions %in% all_exclusions)]
alcohol_exclusions = alcohol_exclusions[!(alcohol_exclusions %in% all_exclusions)]
cvd_exclusions = cvd_exclusions[!(cvd_exclusions %in% all_exclusions)]
stroke_exclusions = stroke_exclusions[!(stroke_exclusions %in% all_exclusions)]
t2d_exclusions = t2d_exclusions[!(t2d_exclusions %in% all_exclusions)]
hypertension_exclusions = hypertension_exclusions[!(hypertension_exclusions %in% all_exclusions)]
glucose_exclusions = glucose_exclusions[!(glucose_exclusions %in% all_exclusions)]
hba1c_exclusions = hba1c_exclusions[!(hba1c_exclusions %in% all_exclusions)]

genText(cancer_exclusions)
genText(hypertension_exclusions)
genText(stroke_exclusions)
genText(cvd_exclusions)
genText(t2d_exclusions)
genText(alcohol_exclusions)
genText(smoking_exclusions)
genText(height_exclusions)
genText(glucose_exclusions)
genText(hba1c_exclusions)

c(numericBinary_key, factors_key)


totalList=c(numericBinary_key$V1, factors_key$V1)
numericsText=""
buffText=""
for (i in 1:length(totalList) ) {
if(i != 1) {buffText=', '}
numericsText = paste0(numericsText, buffText, totalList[i] )
print(i)
}
numericsText



quit()

###########################



# get final exclusion lists
cat $NON_EUR $SEX_DISC_or_LQ $WITHDRAWALS $CLOSE_REATEDS > $rawLoc$'MASTER_EXCLUDE_dupes'
head $rawLoc$'MASTER_EXCLUDE_dupes' 
wc -l  $rawLoc$'MASTER_EXCLUDE_dupes' # 119213


awk '!visited[$0]++' $rawLoc$'MASTER_EXCLUDE_dupes' > $rawLoc$'MASTER_EXCLUDE'
wc -l $rawLoc$'MASTER_EXCLUDE' # 110813

awk 'FNR == NR { lookup[$1] = $1; next; }
FNR <= NR {  
if ($2 in lookup == 0 ) { print $1"\t"$2 } } 
' $rawLoc$'MASTER_EXCLUDE' $UKBB_PLINK1$'ALL.fam' > $rawLoc$'MASTER_KEEP'$Sex_Both
head $rawLoc$'MASTER_KEEP'$Sex_Both
wc -l $rawLoc$'MASTER_KEEP'$Sex_Both # 376761  376,761 indis kept


awk 'FNR == NR { lookup[$1] = $1; next; }
FNR <= NR {  
if ($2 in lookup == 0 && $5 == 1) { print $1"\t"$2 } } 
' $rawLoc$'MASTER_EXCLUDE' $UKBB_PLINK1$'ALL.fam' > $rawLoc$'MASTER_KEEP'$Sex_Male
head $rawLoc$'MASTER_KEEP'$Sex_Male
wc -l $rawLoc$'MASTER_KEEP'$Sex_Male # 174303 Males

awk 'FNR == NR { lookup[$1] = $1; next; }
FNR <= NR {  
if ($2 in lookup == 0 && $5 == 2) { print $1"\t"$2 } } 
' $rawLoc$'MASTER_EXCLUDE' $UKBB_PLINK1$'ALL.fam' > $rawLoc$'MASTER_KEEP'$Sex_Female
head $rawLoc$'MASTER_KEEP'$Sex_Female
wc -l $rawLoc$'MASTER_KEEP'$Sex_Female # 202458 Females

# get sex ratio
#awk '{count[$5]++} END {for (word in count) print word, count[word]}' $UKBB_PLINK1$'ALL.fam'
#1 223033
#2 264362
##############################



# create a PLINK subset of the 375K post QC people, so that we dont have to load in 500K indis into ram in python
arguments=' --bfile '$UKBB_PLINK1$'ALL --memory  43000 --threads 32 --keep '$rawLoc$'MASTER_KEEP'$Sex_Both$'  --keep-allele-order --make-bed --out '$UKBB_PLINK1$'ALL_KEEP'
/home/mk907/software/plink2/plink2 $arguments





################################
j=35
# 4. download PRS files from PGS Catalog and process them ( but only those that we had pheno data for in the UKBB)
# for this I needed to manually upload the 

arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
pheno_binary=${phearray_binary[$j-1]}
# find the sex to be used for pheno
Sex=${sex_subsetList[$j-1]}


rm -rf $rawLoc$pheno_name
mkdir -p $rawLoc$pheno_name
cd $rawLoc$pheno_name

# need to remove windows new line characters...: gsub("\r","",$8)
link=$(awk -v pheno_name="$pheno_name" ' { gsub("\r","",$8); if($8 == pheno_name) {print $4 } }' FS="\t" $rawLoc$'PGSlist_50K_best.txt')
wget $link

# get the filename, by simply listing the dir
filenam=$(ls)
echo $filenam
gunzip $filenam



# get the unzipped filename, by simply listing the dir again
#filenam=$(ls)
# this doesnt work, in case there were multiple files in the zip 

arr=($rawLoc$pheno_name/*) # get the files as an array # https://stackoverflow.com/questions/21668471/bash-script-create-array-of-all-files-in-a-directory
# filenam=${arr[-1]} # get the last element of an array: https://unix.stackexchange.com/questions/198787/is-there-a-way-of-reading-the-last-element-of-an-array-with-bash
filenam=$(basename ${arr[-1]})    # get the filename from path: https://stackoverflow.com/questions/3362920/get-just-the-filename-from-a-path-in-a-bash-script
echo $filenam


# convert it into the format knet expects it

# ignore lines that start with #  https://unix.stackexchange.com/questions/174715/how-to-ignore-the-lines-starts-with-using-grep-awk
#awk -F: '/^[^#]/ {  print $8"\t"$9"\t"$3"\t"$4"\t"$5 }' FS="\t" $filenam > $filenam$knetPRSext
#head $filenam$knetPRSext
#   1       2           3              4               5
# hm_chr  hm_pos  effect_allele   other_allele    effect_weight


# problem is that some PRS files do not have the same columns, but have others... in which case the column numbers wont match:
# 1                  2               3              4                5               6                7                                     8                9            10       11                   12
#rsID            chr_name        chr_position    effect_allele   other_allele    effect_weight   variant_description                     hm_source       hm_rsID         hm_chr  hm_pos       hm_inferOtherAllele
#rs11166389      1               100000723       A               G               -2.302147e-05   FinnGen_VariantID=chr1_100000723_G_A    ENSEMBL         rs11166389      1       100466279

# so we grab the correct columns by name instead

# first get rid of the header to simplify:
awk -F: '/^[^#]/ {  print $0 }' FS="\t" $filenam > $filenam$'_noheader'
# head -n 30 $filenam$'_noheader'

# then perform selection by column name: # https://unix.stackexchange.com/questions/359697/print-columns-in-awk-by-header-name
awk '
NR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
    }
}
{ print $(f["hm_chr"])"\t"$(f["hm_pos"])"\t"$(f["effect_allele"])"\t"$(f["other_allele"])"\t"$(f["effect_weight"]) }' FS="\t" $filenam$'_noheader' > $filenam$knetPRSext
head $filenam$knetPRSext



# remove_ambiguous_alleles: dont use this. As these are not raw GWAS stats, but LD post processed ones, so we could be messing up the LD-adjusted weights 

# create a renamed version of the PGS file for easier loading
cp $rawLoc$pheno_name$'/'$filenam$knetPRSext $rawLoc$pheno_name$'/'$knetPRSext

# also convert it into PLINK PRS format
awk 'FNR == NR { lookup[$1"_"$2] = $3"\t"$5;  forward[$1"_"$2] = $3"_"$4; reverse[$1"_"$2] =  $4"_"$3; next; }
FNR <= NR {  chr_pos=$1"_"$4
if (chr_pos in lookup && ( forward[chr_pos] == $5"_"$6 ||  reverse[chr_pos] == $6"_"$5) ) { print $2"\t"lookup[chr_pos] } } 
' $filenam$knetPRSext $UKBB_PLINK1$'ALL.bim' > $filenam$PLINKPRSext
head $filenam$PLINKPRSext
wc -l $filenam$PLINKPRSext



# create subset of phenotypes to the appropriate sex
awk 'FNR == NR { lookup[$2] = $2; next; } FNR <= NR {  if ($2 in lookup  ) { print $0 } } 
' $rawLoc$'MASTER_KEEP'$Sex $rawLoc$pheno_name$'_all' > $rawLoc$pheno_name$'_all_Sex'
head $rawLoc$pheno_name$'_all_Sex'
wc -l $rawLoc$pheno_name$'_all_Sex'


arguments=' --bfile '$UKBB_PLINK1$'ALL --memory  43000 --threads 32 --score '$rawLoc$pheno_name$'/'$filenam$PLINKPRSext$' ignore-dup-ids cols=-scoreavgs,+scoresums --keep '$rawLoc$'MASTER_KEEP'$Sex$' --pheno '$rawLoc$pheno_name$'_all --out '$rawLoc$pheno_name$'/PRS'
#/home/mk907/software/plink2/plink2 $arguments
sbatch --mem 43000 --cpus-per-task 16 --time 4:0:0 --partition cclake --account danesh-sl3-cpu --wrap "/home/mk907/software/plink2/plink2 $arguments"


#evaluate performance via correlator of the raw PRS
#awk '{print $2"\t"$3"\t"$6 }' $rawLoc$pheno_name$'/PRS.sscore' > $rawLoc$pheno_name$'/PRS_all.sscore'
#Evaluate_PRS $rawLoc$pheno_name$'/PRS_all.sscore' $rawLoc$pheno_name$'/PRS' $pheno_binary '0'#
done # end of pheno extract loop

wc -l $UKBB_PLINK1$'ALL.bim' # 1188672 SNPs in total


#evaluate performance via correlator of the raw PRS
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
pheno_binary=${phearray_binary[$j-1]}

awk '{print $2"\t"$3"\t"$6 }' $rawLoc$pheno_name$'/PRS.sscore' > $rawLoc$pheno_name$'/PRS_all.sscore'
Evaluate_PRS $rawLoc$pheno_name$'/PRS_all.sscore' $rawLoc$pheno_name$'/PRS' $pheno_binary '0' #

done # end of pheno extract loop


#############
# PROBLEM CASES: 
# alcohol consumption, number of cases is virtually everyone
# largest pheno was: EFO_0007878, num cases: 351958
# alcholol-nonconsumption is more interesting, so we swap case/control status
pheno_name="EFO_0007878"
#cp $rawLoc$pheno_name$'_all' $rawLoc$pheno_name$'_all_orig' 
#cp $rawLoc$pheno_name$'_all_Sex' $rawLoc$pheno_name$'_all_Sex_orig'

awk '{if ($3 == 2) {print $1" "$2" 1"} else if ($3 == 1) {print $1" "$2" 2"}  else {print $1" "$2" "$3}}'  $rawLoc$pheno_name$'_all_orig'  >  $rawLoc$pheno_name$'_all'
awk '{if ($3 == 2) {print $1" "$2" 1"} else if ($3 == 1) {print $1" "$2" 2"}  else {print $1" "$2" "$3}}'  $rawLoc$pheno_name$'_all_Sex_orig' >  $rawLoc$pheno_name$'_all_Sex'


#########################

# prioritise controls with complete covariate infos: count how many missing covars are per indi
R
library(readstata13)
library('data.table')
ceu <- read.dta13("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/analysis.dta")
setDT(ceu)
vars = attributes(ceu)$var.labels
names(vars) = names(ceu)

numMissings = rep(0,length(ceu$idno))
for (i in 12:227 ) { # 227# 12 is the age, and 227 is the wage, the rest are non-interesting variables
print(i)
values = ceu[[ colnames(ceu)[i] ]]

if( i >= 106 && i <= 108) {next} # we skip the dates, as these would get converted incorrectly, and we don't need them, as we have these converted to ages as the next few vars
if( i >= 192 && i <= 193) {next}  # also exclude subjecteplist  and subjecteplistf
if( i == 20 || i == 39 || i == 105) {next}  # also exclude married, hxochd, survtype   , as it is all missing

# if it only has missing for more than half then skip
missingVals =  is.na(values) *1
numMissings = numMissings + missingVals

#print(paste0(print(colnames(ceu)[i]), " treatAsNumeric: ", treatAsNumeric)) # vars[i]

}


numMissingsDF=cbind( as.character(ceu$idno),as.numeric( numMissings ) )
png(filename="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/numMissingPerIndi.png"  , width=720, height=720, res=128)
hist(as.numeric(numMissingsDF[,2]))
dev.off()
  
# write to disk
write.table(numMissingsDF, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/numMissingPerIndi", sep="\t", row.names = F, col.names = F, quote = FALSE)

# find a cutoff of max 25 missing values and write those to disk too
numMissingsDF_0 = numMissingsDF[which(as.numeric(numMissingsDF[,2]) < 25),]
nrow(numMissingsDF_0) # 287478, but 30 would give us 415901

write.table(cbind(numMissingsDF_0[,1], numMissingsDF_0[,1]), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/indisWithMaxMissing_25", sep="\t", row.names = F, col.names = F, quote = FALSE)

quit()



##########################

# PROBLEM: 375K indis ~@1mil SNPs takes ~350GB RAM. (and that does not include the peak RAM usage when we need to subset the 350GB and potentially hold 2x as much for a moment)
# SOLUTION: since most phenos have only ~20K case indis, we should just aggregate all the cases + add enough controls to have ~ 4.5 cases / control for each trait
# this should be ~100K in total, which is less than 100GBs of RAM

#PROBLEM2: with so many phenotypes, we cannot use a single PLINK genotype dataset limited to 125K AND keep all the cases as there would be more than 125K....
#=> need multiple genotype subsets... 
#1. for continuous + small #case phenos (< 10K case) use old strategy of combining into a single subset reused for all
#2. large #cases: create a genotype subset for each (capping maximum to be 50% of 125K max, eg for hypertension)



# create a training/valid/test split numbers (these can be reused across all phenos, as these will ALWAYS be 125K, even though the indices will be referring to different indis, as we picked different subsets, but always exactly 125K)
# reproducible random number gen
arguments='/home/mk907/scripts/R/randomNums.R '$maxSamplesize$' 42 '$scratchLoc$'randomNums_split'
Rscript $arguments

# we want a split of 6:2:2 (as the number of cases is quite small for many traits, and we just want to show a difference, not looking for a state of the art performance)
numTrainig=$(awk -v maxSamplesize="$maxSamplesize" 'BEGIN { print int(maxSamplesize * 0.6)}')
numValid=$(awk -v maxSamplesize="$maxSamplesize" 'BEGIN { print int(maxSamplesize * 0.2)}')
numTest=$(awk -v maxSamplesize="$maxSamplesize" 'BEGIN { print int(maxSamplesize * 0.2)}')


awk -v numTrainig="$numTrainig" '{if (FNR <= numTrainig) print $0}' $scratchLoc$'randomNums_split' > $scratchLoc$'randomNums_training'
wc -l $scratchLoc$'randomNums_training' # 75000

awk -v numTrainig="$numTrainig" -v numValid="$numValid" '{if (FNR > numTrainig && FNR <= numTrainig + numValid) print $0}' $scratchLoc$'randomNums_split' > $scratchLoc$'randomNums_valid'
wc -l $scratchLoc$'randomNums_valid' # 25000

awk -v numTrainig="$numTrainig" -v numValid="$numValid" -v numTest="$numTest" '{if (FNR > numTrainig + numValid ) print $0}' $scratchLoc$'randomNums_split' > $scratchLoc$'randomNums_test'
wc -l $scratchLoc$'randomNums_test' # 25000

# create master exclude list for poor quality SNPs
cat $scratchLocAAA$'GENOTYPE_FAILED_QC' $scratchLocAAA$'ALL_FAIL_INFO_09' > $scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC'


###########

# Find out the number of cases for each binary pheno to decide on how to divide them into subsets of phenos
binphenos=""
binphenos_cases=""
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
pheno_binary=${phearray_binary[$j-1]}

#for height we use a different covariates file (the one that excludes height)
if [[ "$pheno_binary" == "1" ]]; then
echo "binary pheno: "$pheno_name
# also keep track of the number of cases for each pheno
binphenos=$binphenos$" "$pheno_name
binphenos_cases=$binphenos_cases$" "$currentCases
fi # end of if binary

done # end of pheno extract loop

echo $binphenos # EFO_0000183 EFO_0000305 EFO_0000537 EFO_0000712 EFO_0000756 EFO_0001663 EFO_0002892 EFO_0003871 EFO_0003877 EFO_0003893 EFO_0004193 EFO_0004286 EFO_0005088 EFO_0005570 EFO_0009259 EFO_0600086 EFO_1000354 EFO_1000657 EFO_1001950 MONDO_0000376 MONDO_0001407 MONDO_0001657 MONDO_0002009 MONDO_0002236 MONDO_0004641 MONDO_0004986 MONDO_0005148 MONDO_0005575 MONDO_0007576 EFO_0007878 EFO_0006527
echo $binphenos_cases # 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560 22560
# then use these in an excel file to manually divide these into the fewest number of subsets, phearray_smallcases, phearray_largecases1, phearray_largecases2, phearray_largecases3




# generate array that holds the information for each phenotype that which PLINK genotype subset it will use for analysis
stringData="phearray_subsets=( "
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}

# by default it will default to the small_cases ( this will work for the quantitative traits too)
pheN='_smallcases'

# if its found on any of the other 
if [[ " ${phearray_largecases1[*]} " =~ " ${pheno_name} " ]]; then
pheN='_largecases1'
fi

if [[ " ${phearray_largecases2[*]} " =~ " ${pheno_name} " ]]; then
pheN='_largecases2'
fi

if [[ " ${phearray_largecases3[*]} " =~ " ${pheno_name} " ]]; then
pheN='_largecases3'
fi

stringData=$stringData$" '"$pheN$"'"

done

stringData=$stringData$" )"

echo $stringData

#######

# create the 4 PLINK subsets
CreatePLINKSubset "_smallcases" "${phearray_smallcases[@]}"

CreatePLINKSubset "_largecases1" "${phearray_largecases1[@]}"

CreatePLINKSubset "_largecases2" "${phearray_largecases2[@]}"

CreatePLINKSubset "_largecases3" "${phearray_largecases3[@]}"


j=1
# need to create phenotype subsets for the train/valid/test splits too, as we cannot rely on getting the appropriate subsets based on the PLINK file, as for scenarios that do NOT use the plink file (eg covar only runs), the valid test would include ALL individuals ( incuding the test sets
# (these phenos are NOT in the same order as the PLINK files, but that is OK, as we match them in python)
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
pheno_binary=${phearray_binary[$j-1]}
pheN=${phearray_subsets[$j-1]}



awk 'FNR == NR { lookup[$2] = $2; next; }
FNR <= NR { if ($2 in lookup  ) { print $0 } } 
' $UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$'.fam' $rawLoc$pheno_name$'_all_Sex' > $rawLoc$pheno_name$'_all_Sex_trainValid'
head $rawLoc$pheno_name$'_all_Sex_trainValid'
wc -l $rawLoc$pheno_name$'_all_Sex_trainValid'

awk 'FNR == NR { lookup[$2] = $2; next; }
FNR <= NR { if ($2 in lookup ) { print $0 } } 
' $UKBB_PLINK1$'ALL_KEEP_TEST'$pheN$'.fam' $rawLoc$pheno_name$'_all_Sex' > $rawLoc$pheno_name$'_all_Sex_test'
head $rawLoc$pheno_name$'_all_Sex_test'
wc -l $rawLoc$pheno_name$'_all_Sex_test'


done




##############################
# II) Analysis (to be run on the GPU nodes)
##############################
$UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN
$UKBB_PLINK1$'ALL_KEEP_TEST'$pheN




##############################
# needed to manually fix the PRS file for EFO_0006525, as it had no chromosome numbers for a few SNPs 
# cp /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS_old

awk '{if ($1 == "T" || $1 == "G" || $1 == "A" || $1 == "C") {print $0}}' /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS_old  > /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS
wc -l /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS # 98 SNPs were bad

awk '{if ($1 == "T" || $1 == "G" || $1 == "A" || $1 == "C") {} else {print $0}}' /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS_old  > /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS
wc -l /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/EFO_0006525/_knetPRS # 1055714
##############################

# (Need to do the below in a separate loop, as the above is too long text for the screen instance to digest)
# E) PRS-only (we just subset the PRS file we generated before to the testset indis
echo -e "pheno\tr2\tp" > $rawLoc$'PRSsigphenos'
j=6
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}

# find the appropriate keep files for PLINK and also the validation sets
pheN=${phearray_subsets[$j-1]}

pheno_binary=${phearray_binary[$j-1]}

# find the appropriate Covars file
covs=${covariateExclusionList[$j-1]}

mkdir -p $resultsLoc$pheno_name$'/eval/'

#awk 'FNR == NR { lookup[$2] = $2; next; } FNR <= NR {  if ($1 in lookup  ) { print $0 } } 
#' $rawLoc$'test_keep'  $rawLoc$pheno_name$'/PRS_all.sscore' > $resultsLoc$pheno_name$'res_PRSindi'
awk 'FNR == NR { lookup[$2] = $2; next; } FNR <= NR {  if ($1 in lookup  ) { print $0 } } 
' $rawLoc$'test_keep'$pheN  $rawLoc$pheno_name$'/PRS_all.sscore' > $resultsLoc$pheno_name$'res_PRSindi'
head $resultsLoc$pheno_name$'res_PRSindi'
wc -l $resultsLoc$pheno_name$'res_PRSindi'

# E) PRS-only (we just subset the PRS file we generated before to the testset indis
awk 'FNR == NR { lookup[$2] = $3; next; } FNR <= NR { if(FNR == 1) {print "ID\tpheno\tPRS"} if ($1 in lookup  ) { print $1"\t"lookup[$1]"\t"$3 } } 
' $rawLoc$pheno_name$'_all_Sex_test' $resultsLoc$pheno_name$'res_PRSindi' >  $resultsLoc$pheno_name$'res_PRSindi_withID2'
head $resultsLoc$pheno_name$'res_PRSindi_withID2'
Evaluate_PRS $resultsLoc$pheno_name$'res_PRSindi_withID2' $resultsLoc$pheno_name$'/eval/res_PRSonly' $pheno_binary '0'

# find out if it was bonf significant, if yes add it to array 
res_PRSindi_p=$(awk '{ if(FNR == 2) {print $3 } }' $resultsLoc$pheno_name$'/eval/res_PRSonly_r2')
#res_PRSindi_p=0.05
isSig=$(awk -v res_PRSindi_p="$res_PRSindi_p" -v arraylength="$arraylength" 'BEGIN { if (res_PRSindi_p < 0.05/arraylength) {print "SIG"} }')

if [ ${#isSig} -ne 0 ] ; then
echo -e $pheno_name$"\t"$res_PRSindi$"\t"$res_PRSindi_p >> $rawLoc$'PRSsigphenos'

# Fit Linear Covar and PRS+Covar models
arguments='/home/mk907/scripts/R/regressionPRSCov.R '$rawLoc$pheno_name$'_all_Sex_trainValid '$rawLoc$'test_keep'$pheN$' '$rawLoc$'valid_keep'$pheN$' '$rawLoc$pheno_name$'/PRS_all.sscore '$COVS_ids$' '$COVS_numericBinary$'_'$covs$' '$COVS_factors$'_'$covs$' '$pheno_binary$' '$resultsLoc$pheno_name$'reg_'
Rscript $arguments
#sbatch --mem 13000 --cpus-per-task 1 --time 2:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

fi

done  # end of submit loop 


# Check if the covariates were all successfully generated
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then

outfile=$resultsLoc$pheno_name$'reg_Cov_CovFactor'
if [ -s "$outfile" ] ; then
bqsd=1
else 
echo "DOES NOT EXIST: "$outfile
fi

fi # if PRS significant

done  # end of submit loop 
# results will be named
#$resultsLoc$pheno_name$'reg_Cov'
#$resultsLoc$pheno_name$'reg_PRS_Cov'


# GENERATE results for the COV and COV+PRSIndi additive baselines
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
pheno_binary=${phearray_binary[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then

#evaluate performance via correlator of the raw PRS
Evaluate_PRS $resultsLoc$pheno_name$'reg_Cov' $resultsLoc$pheno_name$'reg_Cov_res' $pheno_binary '0'
Evaluate_PRS $resultsLoc$pheno_name$'reg_PRS_Cov' $resultsLoc$pheno_name$'reg_PRS_Cov_res' $pheno_binary '0'

fi # if PRS significant
done  # end of submit loop 


# !!! HERE
# results will be named
#$resultsLoc$pheno_name$'reg_Cov_res_r2'
#$resultsLoc$pheno_name$'reg_PRS_Cov_res_r2'



# (this needs to be submitted a few times, due to time limits)
pheno_name="MONDO_0002009"
outName='res_SNP_PRS'
pheno_binary="1"
j=2


# Train KNET:
###############
# training protocol:
# 1. train them with redoEarlyConv=0, keep resubmitting until they all finished (ie not killed via timeout) checking for this via uncommenting res1=$(CheckConvForSubmission2 $outDirr$outName)
# 2. retraing again with redoEarlyConv=1, those that converged within 1 epoch via uncommenting res1=$(CheckConvForSubmission $outDirr$outName)
# 3. to make sure the above finished without killed due to timelimit, apply step (1) again
#######################

redoEarlyConv=0  # when submitting and finishing for the first time
#redoEarlyConv=1 # when all of them finished at least once, and we resubmit to quarter the LR

arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then

pheno_binary=${phearray_binary[$j-1]}
# find the appropriate keep files for PLINK and also the validation sets
pheN=${phearray_subsets[$j-1]}
# find the appropriate Covars file
covs=${covariateExclusionList[$j-1]}

# loop the regular vs _small NN runs
for (( w=0; w<2; w++ )); do
if [ $w == "0" ] ; then 
outDirr=$resultsLoc$pheno_name$'/'
firstLayerSizeUsed=$firstLayerSize
else
outDirr=$resultsLoc$pheno_name$'_small/'
firstLayerSizeUsed=$firstLayerSize_small
fi # end of if first one
echo $outDirr


mkdir -p $outDirr

# A1) SNP+PRS+Covars
outName='res_SNP_PRS_Cov'
arguments='--out  '$outDirr$outName$' knet  --batch_size '$batch_size$' --redoEarlyConv '$redoEarlyConv$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc '$pheno_binary$' --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSizeUsed$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$' --pheno '$rawLoc$pheno_name$'_all_Sex_trainValid  --covars_cont '$resultsLoc$pheno_name$'reg_Cov_CovNumeric --covars_factor '$resultsLoc$pheno_name$'reg_Cov_CovFactor --covars_IDs '$resultsLoc$pheno_name$'reg_Cov_IDs --prs '$rawLoc$pheno_name$'/'$knetPRSext$' --validSet '$rawLoc$'valid_keep'$pheN$' --saveWeights '$outDirr$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
#res1=$(CheckConvForSubmission2 $outDirr$outName)
res1=$(CheckConvForSubmission $outDirr$outName)
if [ ${#res1} -ne 0 ] ; then
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account BUTTERWORTH-SL3-GPU -e $outDirr$outName$'err' -o $outDirr$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
fi


# C1) SNP-only + PRS
outName='res_SNP_PRS'
=arguments='--out  '$outDirr$outName$' knet  --batch_size '$batch_size$'  --redoEarlyConv '$redoEarlyConv$'--gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc '$pheno_binary$' --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSizeUsed$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$' --pheno '$rawLoc$pheno_name$'_all_Sex_trainValid --prs '$rawLoc$pheno_name$'/'$knetPRSext$' --validSet '$rawLoc$'valid_keep'$pheN$' --saveWeights '$outDirr$outName$'SaveWeights' 
#res1=$(CheckConvForSubmission2 $outDirr$outName)
res1=$(CheckConvForSubmission $outDirr$outName)
if [ ${#res1} -ne 0 ] ; then
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account BUTTERWORTH-SL3-GPU -e $outDirr$outName$'err' -o $outDirr$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
fi

# C2) SNP-only 
outName='res_SNP'
arguments='--out  '$outDirr$outName$' knet --disablePRS 1 --prs '$rawLoc$pheno_name$'/'$knetPRSext$'  --batch_size '$batch_size$'  --redoEarlyConv '$redoEarlyConv$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc '$pheno_binary$' --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSizeUsed$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$' --pheno '$rawLoc$pheno_name$'_all_Sex_trainValid --validSet '$rawLoc$'valid_keep'$pheN$' --saveWeights '$outDirr$outName$'SaveWeights' 
#res1=$(CheckConvForSubmission2 $outDirr$outName)
res1=$(CheckConvForSubmission $outDirr$outName)
if [ ${#res1} -ne 0 ] ; then
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account BUTTERWORTH-SL3-GPU -e $outDirr$outName$'err' -o $outDirr$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
fi


done # end of regular and _small loop
fi  # if PRS was significant
done # end of main loop 






# check if all models have successfully run?
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then

# loop the regular vs _small NN runs
for (( w=0; w<2; w++ )); do
if [ $w == "0" ] ; then 
outDirr=$resultsLoc$pheno_name$'/'
firstLayerSizeUsed=$firstLayerSize
else
outDirr=$resultsLoc$pheno_name$'_small/'
firstLayerSizeUsed=$firstLayerSize_small
fi # end of if first one
echo $outDirr

# A1) SNP+PRS+Covars
outName='res_SNP_PRS_Cov'
CheckIfNonLin_And_Lin $outDirr$outName$'SaveWeights'
CheckIfConvergedLin_And_Lin $outDirr$outName

# C1) SNP-only + PRS
outName='res_SNP_PRS'
#CheckIfNonLin_And_Lin $outDirr$outName$'SaveWeights'
CheckIfConvergedLin_And_Lin $outDirr$outName

# C2) SNP-only 
outName='res_SNP'
#CheckIfNonLin_And_Lin $outDirr$outName$'SaveWeights'
CheckIfConvergedLin_And_Lin $outDirr$outName


done # end of regular and _small loop

fi  # if PRS was significant
done # end of main loop 


############################################

j=1
# Inference KNET:
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then

pheno_binary=${phearray_binary[$j-1]}
# find the appropriate PLINK file
pheN=${phearray_subsets[$j-1]}
# find the appropriate Covars file
covs=${covariateExclusionList[$j-1]}

# loop the regular vs _small NN runs
for (( w=0; w<2; w++ )); do
if [ $w == "0" ] ; then 
outDirr=$resultsLoc$pheno_name$'/'
firstLayerSizeUsed=$firstLayerSize
else
outDirr=$resultsLoc$pheno_name$'_small/'
firstLayerSizeUsed=$firstLayerSize_small
fi # end of if first one
echo $outDirr

# A1) SNP+PRS+Covars
outName='res_SNP_PRS_Cov'
arguments='--out  '$outDirr$outName$' knet --inference 1 --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc '$pheno_binary$' --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSizeUsed$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'ALL_KEEP_TEST'$pheN$' --pheno '$rawLoc$pheno_name$'_all_Sex_test  --covars_cont '$resultsLoc$pheno_name$'reg_Cov_CovNumeric --covars_factor '$resultsLoc$pheno_name$'reg_Cov_CovFactor --covars_IDs '$resultsLoc$pheno_name$'reg_Cov_IDs --prs '$rawLoc$pheno_name$'/'$knetPRSext$' --loadWeights '$outDirr$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 2:0:0 --partition ampere --account BUTTERWORTH-SL3-GPU -e $outDirr$outName$'_inf_err' -o $outDirr$outName$'_inf_out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"


# C1) SNP-only + PRS
outName='res_SNP_PRS'

arguments='--out  '$outDirr$outName$' knet  --inference 1 --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc '$pheno_binary$' --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSizeUsed$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'ALL_KEEP_TEST'$pheN$' --pheno '$rawLoc$pheno_name$'_all_Sex_test --prs '$rawLoc$pheno_name$'/'$knetPRSext$' --loadWeights '$outDirr$outName$'SaveWeights' 
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 2:0:0 --partition ampere --account BUTTERWORTH-SL3-GPU -e $outDirr$outName$'_inf_err' -o $outDirr$outName$'_inf_out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"

# C2) SNP-only 
outName='res_SNP'
arguments='--out  '$outDirr$outName$' knet  --inference 1 --disablePRS 1 --prs '$rawLoc$pheno_name$'/'$knetPRSext$'  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc '$pheno_binary$' --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSizeUsed$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'ALL_KEEP_TEST'$pheN$' --pheno '$rawLoc$pheno_name$'_all_Sex_test --loadWeights '$outDirr$outName$'SaveWeights' 
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 2:0:0 --partition ampere --account BUTTERWORTH-SL3-GPU -e $outDirr$outName$'_inf_err' -o $outDirr$outName$'_inf_out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"

# E) PRS-only (we just subset the PRS file we generated before to the testset indis
done # end of regular and _small loop
fi  # if PRS was significant
done # end of main loop 



# $UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$'.bim'
# $UKBB_PLINK1$'ALL_KEEP_TEST'$pheN$'.bim'

# # the output files for the above
# $outDirr$outName$"FIDs_TEST.txt"
# $outDirr$outName$"yhat_TEST_noAct_retrain.txt"
# $outDirr$outName$"yhat_TEST_noAct.txt"
# $outDirr$outName$"yhat_TEST.txt"


# !!! HERE
# Check if all inference succesfully run
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then


# loop the regular vs _small NN runs
for (( w=0; w<2; w++ )); do
if [ $w == "0" ] ; then 
outDirr=$resultsLoc$pheno_name$'/'
else
outDirr=$resultsLoc$pheno_name$'_small/'
fi # end of if first one
echo $outDirr

# A1) SNP+PRS+Covars
outName='res_SNP_PRS_Cov'
Check_TEST_exists $outDirr$outName

# C1) SNP-only + PRS
outName='res_SNP_PRS'
Check_TEST_exists $outDirr$outName

# C2) SNP-only 
outName='res_SNP'
Check_TEST_exists $outDirr$outName


done # end of regular and _small loop
fi  # if PRS was significant
done # end of main loop 





# Evaluate PRS
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then

pheno_binary=${phearray_binary[$j-1]}

# loop the regular vs _small NN runs
for (( w=0; w<2; w++ )); do
if [ $w == "0" ] ; then 
outDirr=$resultsLoc$pheno_name$'/'
firstLayerSizeUsed=$firstLayerSize
else
outDirr=$resultsLoc$pheno_name$'_small/'
firstLayerSizeUsed=$firstLayerSize_small
fi # end of if first one
echo $outDirr

mkdir -p $outDirr$'eval/'


# A1) SNP+PRS+Covars
outName='res_SNP_PRS_Cov'
GenAndEval_ALL_PRS $outDirr$outName $rawLoc$pheno_name$'_all_Sex_test' $pheno_binary $outDirr$'eval/'$outName '0'

# C1) SNP-only + PRS
outName='res_SNP_PRS'
GenAndEval_ALL_PRS $outDirr$outName $rawLoc$pheno_name$'_all_Sex_test' $pheno_binary $outDirr$'eval/'$outName '0'

# C2) SNP-only 
outName='res_SNP'
GenAndEval_ALL_PRS $outDirr$outName $rawLoc$pheno_name$'_all_Sex_test' $pheno_binary $outDirr$'eval/'$outName '0'

done # end of regular and _small loop
fi  # if PRS was significant
done # end of main loop 





# BUILD FINAL RESULTS TABLE
# loop the regular vs _small NN runs
for (( w=0; w<2; w++ )); do
if [ $w == "0" ] ; then 
smallPostfix=''
else
smallPostfix='_small'
fi # end of if first one

# extract all results into 3 text files ( we work off the correlation not correlation^2)
header="pheno\tSNP+PRS+Cov\tSNP+PRS\tSNP\tPRSindi\tPRSindi+CovAdd"
echo -e $header > $resultsLoc$'finalRes_nonlinear'$smallPostfix
echo -e $header > $resultsLoc$'finalRes_noAct'$smallPostfix
echo -e $header > $resultsLoc$'finalRes_linear'$smallPostfix

#header="pheno\tPRS_p"
echo -e $header > $resultsLoc$'finalRes_nonlinear_p'$smallPostfix
echo -e $header > $resultsLoc$'finalRes_noAct_p'$smallPostfix
echo -e $header > $resultsLoc$'finalRes_linear_p'$smallPostfix

# generate results for the Validation set too:
header="pheno\tSNP+PRS+Cov\tSNP+PRS\tSNP"
echo -e $header > $resultsLoc$'finalRes_nonlinear_valid'$smallPostfix
echo -e $header > $resultsLoc$'finalRes_linear_valid'$smallPostfix
done # end of regular and _small loop


j=8
arraylength=${#phearray_withCEUdef[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${phearray_withCEUdef[$j-1]}
onTheList=$(grep $pheno_name $rawLoc$'PRSsigphenos')
if [ ${#onTheList} -ne 0 ] ; then

# loop the regular vs _small NN runs
for (( w=0; w<2; w++ )); do
if [ $w == "0" ] ; then 
smallPostfix=''
outDirr=$resultsLoc$pheno_name$'/'
else
smallPostfix='_small'
outDirr=$resultsLoc$pheno_name$'_small/'
fi # end of if first one

echo $outDirr

mkdir -p $outDirr$'eval/'

# A1) SNP+PRS+Covars
outName='res_SNP_PRS_Cov'
res_SNP_PRS_Cov_linear=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_linear_r2')
res_SNP_PRS_Cov_linear_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_linear_r2')
res_SNP_PRS_Cov_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_nonlinear_r2')
res_SNP_PRS_Cov_nonlinear_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_nonlinear_r2')
res_SNP_PRS_Cov_noAct=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_noAct_r2')
res_SNP_PRS_Cov_noAct_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_noAct_r2')
res_SNP_PRS_Cov_nonlinear_valid=$(awk '{ print $1  }' $outDirr$outName$'_bestKNET_PRS')
res_SNP_PRS_Cov_linear_valid=$(awk '{ print $1  }' $outDirr$outName$'_linear_bestKNET_PRS')


# C1) SNP-only + PRS
outName='res_SNP_PRS'
res_SNP_PRS_linear=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_linear_r2')
res_SNP_PRS_linear_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_linear_r2')
res_SNP_PRS_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_nonlinear_r2')
res_SNP_PRS_nonlinear_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_nonlinear_r2')
res_SNP_PRS_noAct=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_noAct_r2')
res_SNP_PRS_noAct_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_noAct_r2')
res_SNP_PRS_nonlinear_valid=$(awk '{ print $1  }' $outDirr$outName$'_bestKNET_PRS')
res_SNP_PRS_linear_valid=$(awk '{ print $1  }' $outDirr$outName$'_linear_bestKNET_PRS')

# C2) SNP-only 
outName='res_SNP'
res_SNP_linear=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_linear_r2')
res_SNP_linear_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_linear_r2')
res_SNP_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_nonlinear_r2')
res_SNP_nonlinear_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_nonlinear_r2')
res_SNP_noAct=$(awk '{ if(FNR == 2) {print $1 } }'  $outDirr$'eval/'$outName$'_noAct_r2')
res_SNP_noAct_p=$(awk '{ if(FNR == 2) {print $3 } }'  $outDirr$'eval/'$outName$'_noAct_r2')
res_SNP_nonlinear_valid=$(awk '{ print $1  }' $outDirr$outName$'_bestKNET_PRS')
res_SNP_linear_valid=$(awk '{ print $1  }' $outDirr$outName$'_linear_bestKNET_PRS')


# E1) PRS-only (we just subset the PRS file we generated before to the testset indis
res_PRSindi=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$pheno_name$'/eval/res_PRSonly_r2')
res_PRSindi_p=$(awk '{ if(FNR == 2) {print $3 } }' $resultsLoc$pheno_name$'/eval/res_PRSonly_r2')

# E2) PRS + Cov additive baseline
res_PRSindi_CovAdd=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$pheno_name$'reg_PRS_Cov_res_r2')
res_PRSindi_CovAdd_p=$(awk '{ if(FNR == 2) {print $3 } }' $resultsLoc$pheno_name$'reg_PRS_Cov_res_r2')


echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_nonlinear$"\t"$res_SNP_PRS_nonlinear$"\t"$res_SNP_nonlinear$"\t"$res_PRSindi$"\t"$res_PRSindi_CovAdd >> $resultsLoc$'finalRes_nonlinear'$smallPostfix
echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_linear$"\t"$res_SNP_PRS_linear$"\t"$res_SNP_linear$"\t"$res_PRSindi$"\t"$res_PRSindi_CovAdd >> $resultsLoc$'finalRes_linear'$smallPostfix
echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_noAct$"\t"$res_SNP_PRS_noAct$"\t"$res_SNP_noAct$"\t"$res_PRSindi$"\t"$res_PRSindi_CovAdd >> $resultsLoc$'finalRes_noAct'$smallPostfix

# also get the p-val for the PRS
echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_nonlinear_p$"\t"$res_SNP_PRS_nonlinear_p$"\t"$res_SNP_nonlinear_p$"\t"$res_PRSindi_p$"\t"$res_PRSindi_CovAdd_p >> $resultsLoc$'finalRes_nonlinear_p'$smallPostfix
echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_linear_p$"\t"$res_SNP_PRS_linear_p$"\t"$res_SNP_linear_p$"\t"$res_PRSindi_p$"\t"$res_PRSindi_CovAdd_p >> $resultsLoc$'finalRes_linear_p'$smallPostfix
echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_noAct_p$"\t"$res_SNP_PRS_noAct_p$"\t"$res_SNP_noAct_p$"\t"$res_PRSindi_p$"\t"$res_PRSindi_CovAdd_p >> $resultsLoc$'finalRes_noAct_p'$smallPostfix

# generate results for the Validation set too:
echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_nonlinear_valid$"\t"$res_SNP_PRS_nonlinear_valid$"\t"$res_SNP_nonlinear_valid >> $resultsLoc$'finalRes_nonlinear_valid'$smallPostfix
echo -e $pheno_name$"\t"$res_SNP_PRS_Cov_linear_valid$"\t"$res_SNP_PRS_linear_valid$"\t"$res_SNP_linear_valid >> $resultsLoc$'finalRes_linear_valid'$smallPostfix


done # end of regular and _small loop

fi  # if PRS was significant
done # end of main loop 




#############################################
# VI. Plot Results
#############################################


arguments='/home/mk907/scripts/R/PRS_GXE_Plotter.R '$resultsLoc$'finalRes '$resultsLoc$'finalRes_res PRS_GxE_non-linearity'
Rscript $arguments

# "SNP linear to non-linear % difference: 0.0686 of 27 traits"
# "SNP+Cov linear to non-linear % difference: 0.0377 of 21 traits"
# "SNP non-linear to baselin % difference: 0.0704"
# "SNP+Cov non-linear to baseline % difference: 0.0525"

# So it seems that there is again, some, minimal amount of non-linearity, but it is NOT better than the basic, additive regression model

arguments='/home/mk907/scripts/R/PRS_GXE_Plotter.R '$resultsLoc$'0RESULTS/finalRes '$resultsLoc$'finalRes_res_chris PRS_GxE_non-linearity'
/usr/local/Cluster-Apps/R/4.3.1-icelake/bin/Rscript $arguments



#####################


# TODO?:

# Run simulations of the epistatic scenario with the PRS SNP weighting strategy enabled
# this would be needed, as we have seen that even in the case of complete coverage and no genuine epistasis, the nonlinear method was better by 3.6%?
# ( so it seems likely that 11%, is too high
# in the meanwhile remove from the table this result, as this is just confusing

# actually I've done this: res_epi_PRS