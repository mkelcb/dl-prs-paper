############################# 

# screen -r -D 124494.pts-46.login-p-1

# Sims:
# screen -r -D 167694.pts-0.login-p-2

# GPU
# real data:
# screen -r -D 10750.pts-51.login-e-2

# Sims
# screen -r -D 5125.pts-16.login-e-4

# scontrol show jobid -dd 36062877
# reusable functions and variables used by other scripts
#####################################



ukbb_raw_pheno='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/output/AAD/AAD_events_and_followup.txt'
homeScratch='/rds/user/mk907/hpc-work/scratch/'
endpointLoc='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/'
baseLoc='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/'
dataLoc=$baseLoc$'data/'
rawLoc=$dataLoc$'raw/'
scratchLoc=$baseLoc$'scratch/'
resultsLoc=$baseLoc$'results/'
plink="/home/mk907/software/plink/plink"
plink2='/home/mk907/software/plink2_new/plink2'  #'/home/mk907/software/plink2/plink2'
shaPRS="_shaPRS"
NCORES_LDPred2=16
QCfile='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/reference_files/ukb_sqc_v2.txt' # eid matched to our app
# signature:
#  1              2                 3          4            5         6            7               8                   9                        10                  11              12                 13                     14                 15                  16                    17                       18                              19                              20                                           21          22                              23                                24 
# eid     genotyping.array        Batch   Plate.Name      Well    Cluster.CR      dQC     Internal.Pico..ng.uL.   Submitted.Gender        Inferred.Gender       X.intensity     Y.intensity     Submitted.Plate.Name    Submitted.Well  sample.qc.missing.rate  heterozygosity  heterozygosity.pc.corrected  het.missing.outliers     putative.sex.chromosome.aneuploidy      in.kinship.table        excluded.from.kinship.inference excess.relatives        in.white.British.ancestry.subset      used.in.pca.calculation PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10    PC11    PC12 PC13     PC14    PC15    PC16    PC17    PC18    PC19    PC20    PC21    PC22    PC23    PC24    PC25    PC26    PC27    PC28    PC29    PC30    PC31 PC32     PC33    PC34    PC35    PC36    PC37    PC38    PC39    PC40    in.Phasing.Input.chr1_22        in.Phasing.Input.chrX   in.Phasing.Input.chrXY


# common Keeplists shared with other projects
baseLoc_common='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/'
dataLoc_common=$baseLoc_common$'data/'
rawLoc_common=$dataLoc_common$'raw/'
NON_EUR=$rawLoc_common$'NON_EUR' # those who are NOT European ancestry
SEX_DISC_or_LQ=$rawLoc_common$'SEX_DISC_or_LQ'  # Sex discordant or those with low quality genotypes ( too many missing or excess heterozygosity)
META_GRS=$rawLoc_common$'META_GRS' # the 12K training set for the Inouye lab's meta GRS project
WITHDRAWALS=$rawLoc_common$'WITHDRAWALS' # people no longer in the UKBB
#TEST_SET=$rawLoc_common$'EUR_subset' # those that are in the test set (the UKBiLEVE chip)
RELATEDS=$rawLoc_common$'RELATEDS' # pairs of people who are too closely related 
CLOSE_REATEDS=$rawLoc_common$'CLOSE_REATEDS' # the individual to exclude from the above (preferentially controls)
hapmap3_b37bim=$dataLoc_common$'raw/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim'
hapmap3_SNPs=$dataLoc_common$'raw/hapmap3_SNPs'
UKBB_PLINK1=$dataLoc_common$'UKBB_PLINK1/'
endpointLoc='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/'
knetPRSext="_knetPRS"
PLINKPRSext="_PLINKPRS"
scratchLocAAA='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/scratch/'

Sex_Both="all"
Sex_Male="male"
Sex_Female="female"


# how many individuals we could have in an analysis, based on RAM availibility on the nodes
maxSamplesize=125000
maxAllowedCases=62500 # this is half of the above


# subset lists
phearray=( 'EFO_0000183' 'EFO_0000305' 'EFO_0000537' 'EFO_0000712' 'EFO_0000756' 'EFO_0001663' 'EFO_0002892' 'EFO_0003871' 'EFO_0003877' 'EFO_0003893' 'EFO_0004193' 'EFO_0004286' 'EFO_0005088' 'EFO_0005570' 'EFO_0009259' 'EFO_0600086' 'EFO_1000354' 'EFO_1000657' 'EFO_1001950' 'MONDO_0000376' 'MONDO_0001407' 'MONDO_0001657' 'MONDO_0002009' 'MONDO_0002236' 'MONDO_0004641' 'MONDO_0004986' 'MONDO_0005148' 'MONDO_0005575' 'MONDO_0007576' )
phearray_noT2D=( 'EFO_0000183' 'EFO_0000305' 'EFO_0000537' 'EFO_0000712' 'EFO_0000756' 'EFO_0001663' 'EFO_0002892' 'EFO_0003871' 'EFO_0003877' 'EFO_0003893' 'EFO_0004193' 'EFO_0004286' 'EFO_0005088' 'EFO_0005570' 'EFO_0009259' 'EFO_0600086' 'EFO_1000354' 'EFO_1000657' 'EFO_1001950' 'MONDO_0000376' 'MONDO_0001407' 'MONDO_0001657' 'MONDO_0002009' 'MONDO_0002236' 'MONDO_0004641' 'MONDO_0004986' 'MONDO_0005575' 'MONDO_0007576' )

# full pheno lists, with indices aligned
phearray_withCEUdef=( 'EFO_0000183' 'EFO_0000305' 'EFO_0000537' 'EFO_0000712' 'EFO_0000756' 'EFO_0001663' 'EFO_0002892' 'EFO_0003871' 'EFO_0003877' 'EFO_0003893' 'EFO_0004193' 'EFO_0004286' 'EFO_0005088' 'EFO_0005570' 'EFO_0009259' 'EFO_0600086' 'EFO_1000354' 'EFO_1000657' 'EFO_1001950' 'MONDO_0000376' 'MONDO_0001407' 'MONDO_0001657' 'MONDO_0002009' 'MONDO_0002236' 'MONDO_0004641' 'MONDO_0004986' 'MONDO_0005148' 'MONDO_0005575' 'MONDO_0007576' 'EFO_0007878' 'EFO_0006527' 'EFO_0006525' 'EFO_0004339' 'EFO_0004465' 'EFO_0004541' )
phearray_binary=( '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '0' '0' '0' '0' )
sex_subsetList=( 'all' 'female' 'all' 'all' 'all' 'male' 'all' 'all' 'all' 'all' 'all' 'all' 'male' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' 'all' )
covariateExclusionList=( 'cancer' 'cancer' 'hypertension' 'stroke' 'cancer' 'cancer' 'cancer' 'cancer' 'all' 'cancer' 'cancer' 'cvd' 'cancer' 'cancer' 'cancer' 'all' 'cancer' 'cancer' 'cancer' 'cancer' 'cancer' 'cancer' 'all' 'cancer' 'cancer' 'cancer' 'diabetes' 'cancer' 'cancer' 'alcohol' 'smoking' 'smoking' 'height' 'glucose' 'hba1c' )


# 
phearray_smallcases=( 'EFO_0000183' 'EFO_0000756' 'EFO_0002892' 'EFO_0003871' 'EFO_0003877' 'EFO_0003893' 'EFO_0004286' 'EFO_0005088' 'EFO_0005570' 'EFO_0009259' 'EFO_0600086' 'EFO_1000354' 'EFO_1000657' 'EFO_1001950' 'MONDO_0000376' 'MONDO_0001407' 'MONDO_0001657' 'MONDO_0002236' 'MONDO_0004641' 'MONDO_0004986' 'MONDO_0005575' 'MONDO_0007576' )
phearray_largecases1=( 'EFO_0000305' 'EFO_0000712' 'EFO_0001663' 'EFO_0004193' )
phearray_largecases2=( 'EFO_0006527' 'MONDO_0002009' )
phearray_largecases3=( 'EFO_0000537' )


phearray_subsets=( '_smallcases' '_largecases1' '_largecases3' '_largecases1' '_smallcases' '_largecases1' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_largecases1' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_largecases2' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_smallcases' '_largecases2' '_smallcases' '_smallcases' '_smallcases' '_smallcases' )


# sims
simsDataTruthset=$rawLoc$'truthset/'
recombinationBlokcs='/home/mk907/scripts/R/all.1cM.tab'
simsLoc=$rawLoc$'sims/'
simsData=$simsLoc$'data/'

numSNPsToPick=500000
causalPerc="0.004"
h2level="0.5"
numSimReplicates=20
numOrders=4

# 500000 * 0.004 =  2000
halfCausals=1000 # as bash doesnt support floating point properly
echo $halfCausals



j=12

#wc -l $UKBB_PLINK1$'ALL.bim' # 1,188,672
#wc -l $UKBB_PLINK1$'ALL.fam' # 487395
# so this is the full UKB, with all people on hm3...

################
# take a 1000 people, with 10K


#################


numChroms=22
plinkformat='.phe.glm.logistic.hybrid'
plinkformatlinear=".phe.glm.linear"

GENOTYPE_QC='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/pre_qc_data/affy_ukbiobank_array/raw_data/showcase_release_19Jul/snpqc/ukb_snp_qc.txt'

shaPRSscriptLoc="/home/mk907/scripts/R/shaPRS.R"



ldscDir='/home/mk907/software/ldsc/'
w_hm3_snplist='/home/mk907/software/ldsc/w_hm3.snplist'
eur_w_ld_chr='/home/mk907/software/ldsc/eur_w_ld_chr/'
lsdc=$ldscDir$'ldsc/ldsc.py'
munge_sumstats=$ldscDir$'ldsc/munge_sumstats.py'

numChroms=22

PLINKRAM=256000
SLURM_CPUS_ON_NODE=32










mtagscript='/home/mk907/software/mtag/'











# check for rG 
module load miniconda/2
source activate ldsc
module load R/4.0.3
# pip install pandas --user
# pip install scipy --user
# pip install bitarray --user







##############################################
# GPU


COVS_ids="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_ids"
COVS_numericBinary="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary"
COVS_numericBinary_key="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary_key"


COVS_factors="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors"
COVS_factors_key="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors_key"
hidAct='4'
hidCount='3'
firstLayerSize='100'
firstLayerSize_small='24'
dropout='0.3'
gpu='1'
batch_size='32'
gradient_batch_size=$batch_size



#pip3 install hyperopt
#pip3 install prettytable

#pip3 install PyPlink

# pip3 install statsmodels

module load python/3.6 cuda/11.2 cudnn/8.1_cuda-11.2


numChroms=22
module load R/4.0.3

#source activate martin
source ~/pytorch-env113/bin/activate


#python3
#import torch
#print(torch.__version__) # 1.9.0+cu111





#########################################

# exports single pheno from multiple codes, also filters out duplicates
function export_UKBB_pheno_nodupes3 { 
pheno_name=$1
definition=$2
rawLoc=$3

# load the right version of R that has all the packages
module load R/4.0.3

# extraction script only works from this location
endpointLoc='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/'

#endpointLoc='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/deprecated/endpoints/'

# must be in this folder, otherwise script wont work...
cd $endpointLoc
outLoc=$endpointLoc$'output/'$pheno_name
mkdir -p $outLoc

# run script to extract pheno
Rscript /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/curate_endpoint.R  --def-file $definition --output $outLoc

ukbb_raw_pheno=$outLoc'/events_and_followup.txt'
# signature:
#  1        2         3                  4                     5                    6                  7                  8                    9                          10                      11                     12                           13                                   14                             15                            16                    17                        18                     19                      20                      21                      22                         23                   24                         25                             26                           27                      28                       29                   30                   31                                32                          33                              34                               35                       36                              37                    38                   39                40
#eid     sex     visit_index     assessment_date      assessment_centre       assessment_nation       age       ehr_linkage_withdrawn   earliest_hospital_date  earliest_hospital_nation  prevalent_event    prevalent_event_with_followup   prevalent_event_followup        prevalent_event_followup_date   prevalent_event_followup_age      prevalent_event_date    prevalent_event_age     prevalent_event_type    prevalent_record_source prevalent_cause_type    prevalent_code_type       prevalent_code        prevalent_code_label    lost_to_followup        lost_to_followup_reason          lost_to_followup_date       latest_hospital_date    latest_hospital_nation    latest_mortality_date   incident_event  incident_event_followup        incident_event_followup_date    incident_event_followup_age     mortality_at_followup_date        incident_event_type     incident_record_source          incident_cause_type     incident_code_type      incident_code   incident_code_label
#1000012 Female         0           2010-03-17            Liverpool              England              58.25           FALSE                   1993-07-27                    England            FALSE                  FALSE                     -16.6408281998631                   1993-07-27                           41.61                        FALSE                 2021-09-30           England 2021-10-15      FALSE   11.5355920602327  2021-09-30      69.79   FALSE

#awk 'FNR <= NR {  if( FNR < 3 ) {print $40 } }' $ukbb_raw_pheno 
#awk 'FNR <= NR {  if( FNR < 20 ) {print $11"\t"$30 } }' FS="\t" $ukbb_raw_pheno 

# filter duplicates, we want the baseline entry for each individual only IE, where visit == 0, (otherwise the inci_followup will be from their current visit which is not the total)

awk '
FNR <= NR {  
if( FNR == 1 || $3 == 0 ) {print $0 } 
}
' $ukbb_raw_pheno > $rawLoc$pheno_name$'_survival'

original=$(wc -l $ukbb_raw_pheno)
noDupes=$(wc -l $rawLoc$pheno_name$'_survival')
#echo "before/after filtering for dupes: "$original$" / "$noDupes

outExtract=$endpointLoc$'output/'$pheno_name'/'$pheno_name$'_filtered'
#awk '{ if ($11 == "TRUE" || $30 == "TRUE" ) { print $1"\t"$1"\t2" } else {print $1"\t"$1"\t1"}  }' FS="\t" $rawLoc$pheno_name$'_survival' > $rawLoc$pheno_name$'_all'
awk '{ if(FNR > 1)  { if ($11 == "TRUE" || $30 == "TRUE" ) { print $1"\t"$1"\t2" } else {print $1"\t"$1"\t1"} } }' FS="\t" $rawLoc$pheno_name$'_survival' > $rawLoc$pheno_name$'_all'


awk '{count[$3]++} END {for (word in count) print word, count[word]}' $rawLoc$pheno_name$'_all'
}
export -f export_UKBB_pheno_nodupes3 # this makes local functions executable when bsubbed




#allScoreLoc=$rawLoc$pheno_name$'/PRS.sscore'
#outResultLoc=$rawLoc$pheno_name$'/PRS'
#binary=1
# evaluates the r^2 and AUC of a PRS on a test set
function Evaluate_PRS { 
allScoreLoc=$1
outResultLoc=$2
binary=$3
append=$4

arguments='/home/mk907/scripts/R/correlator_r2redux.R '$allScoreLoc$' '$outResultLoc$'_r2 '$binary$' 1 '$append$' 1'
Rscript $arguments
}
export -f Evaluate_PRS # this makes local functions executable when bsubbed


#######################################################

#outFolder=$resultsLoc$pheno_name$'/eval/'$outName
# generates all 3 types of DL models (nonlinear, no-act and retrained no-act)
function GenAndEval_ALL_PRS { 
fil=$1
testpheno=$2
pheno_binary=$3
outFolder=$4
append=$5

GenAndEval_PRS $fil$"FIDs_TEST.txt" $fil$"yhat_TEST_noAct_retrain.txt" $testpheno $pheno_binary $outFolder$'_linear' $append

GenAndEval_PRS $fil$"FIDs_TEST.txt" $fil$"yhat_TEST_noAct.txt" $testpheno $pheno_binary $outFolder$'_noAct' $append

GenAndEval_PRS $fil$"FIDs_TEST.txt" $fil$"yhat_TEST.txt" $testpheno $pheno_binary $outFolder$'_nonlinear' $append

}


# inFile=$resultsLoc$pheno_name$'/eval/res_PRSonly_r2'
function getCor { 
inFile=$1
corr=$(awk '{ if(FNR == 2) {print $1 } }' $inFile)
echo $corr
}



#FIDFile=$fil$"FIDs_TEST.txt"
#yhatFile=$resultsLoc$pheno_name$'/'$outName$"yhat_TEST.txt"
#testpheno=$rawLoc$pheno_name$'_all_Sex_test'
#binary=1
#outF=$outFolder$'_nonlinear'
function GenAndEval_PRS { 
FIDFile=$1
yhatFile=$2
testpheno=$3
binary=$4
outF=$5
append=$6

#paste together the FID and yhatfiles
paste $FIDFile $yhatFile > $yhatFile$'_withID'

# createa an input file format suitable for Evaluate_PRS
awk 'FNR == NR { lookup[$2] = $3; next; } FNR <= NR { if(FNR == 1) {print "ID\tpheno\tPRS"} if ($2 in lookup  ) { print $2"\t"lookup[$2]"\t"$3 } } 
' $testpheno  $yhatFile$'_withID' > $yhatFile$'_withID2'

Evaluate_PRS $yhatFile$'_withID2' $outF $binary $append

}


function CheckConvForSubmission { 
target=$1
res1=$(CheckIfNonLin_And_Lin $target$'SaveWeights')
res2=$(CheckIfConvergedLin_And_Lin $target)
if [ ${#res1} -eq 0 ] && [ ${#res2} -eq 0 ] ; then
#echo "CONVERGED WELL"
uba=1
else
echo "DID NOT CONVERGE"
fi
}

# same as above, but it only checks if its finished training, NOT if its converged within 1 epoch
function CheckConvForSubmission2 { 
target=$1
res1=$(CheckIfNonLin_And_Lin $target$'SaveWeights')
if [ ${#res1} -eq 0 ] ; then
#echo "CONVERGED WELL"
uba=1
else
echo "DID NOT FINISH"
fi
}



# check if the model actually converged (ie best epoch > 1)
# fil=$resultsLoc$pheno_name$'/'$outName
function CheckIfConvergedLin_And_Lin { 
fil=$1
CheckIfConverged $fil$'best_epoch.txt'
CheckIfConverged $fil$'_linearbest_epoch.txt'
}

# outfile=$resultsLoc$pheno_name$'/'$outName$'best_epoch.txt'
# outfile=$resultsLoc$pheno_name$'/'$outName$'_linearbest_epoch.txt'
function CheckIfConverged { 
outfile=$1
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
source $outfile
if (( $best_epoch <= 0 )); then
    echo "DID NOT CONVERGE ("$best_epoch$"): "$outfile
fi
fi
}




function CheckIfNonLin_And_Lin { 
fil=$1

CheckIfExsists $fil

CheckIfExsists $fil$'_linear'
}


function CheckIfExsists { 
outfile=$1
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
bqq=1
else
echo "DOES NOT EXIST: "$outfile
fi
}


function Check_TEST_exists { 
fil=$1

CheckIfExsists $fil$"FIDs_TEST.txt"

CheckIfExsists $fil$"yhat_TEST_noAct_retrain.txt"

CheckIfExsists $fil$"yhat_TEST_noAct.txt"

CheckIfExsists $fil$"yhat_TEST.txt"

}





function DELETE_TEST { 
fil=$1

rm -rf $fil$"FIDs_TEST.txt"

rm -rf $fil$"yhat_TEST_noAct_retrain.txt"

rm -rf $fil$"yhat_TEST_noAct.txt"

rm -rf $fil$"yhat_TEST.txt"

}



# https://askubuntu.com/questions/674333/how-to-pass-an-array-as-function-argument
# loop through and creates subset of PLINK genotype data that includes all the cases and the correct number of controls with least amount of missing data
function CreatePLINKSubset() {
pheN="$1"   # Save first argument in a variable
shift            # Shift all arguments to the left (original $1 gets lost)
pheArray_current=("$@") # Rebuild the array with rest of arguments

# go through all binary traits
rm -rf $scratchLoc$'all_cases_dupes'$pheN
arraylength=${#pheArray_current[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
pheno_name=${pheArray_current[$j-1]}

echo "binary pheno: "$pheno_name

# wc -l $rawLoc$pheno_name$'_all_Sex'
# get all the cases
awk '{if ($3 == 2) {print $1"\t"$2}}' $rawLoc$pheno_name$'_all_Sex' > $scratchLoc$pheno_name$'_cases'
wc -l $scratchLoc$pheno_name$'_cases'

cat $scratchLoc$pheno_name$'_cases' >> $scratchLoc$'all_cases_dupes'$pheN

done # end of pheno extract loop

# remove duplicates 
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $scratchLoc$'all_cases_dupes'$pheN > $scratchLoc$'all_cases'$pheN
head $scratchLoc$'all_cases'$pheN
wc -l $scratchLoc$'all_cases'$pheN # 48798
wc -l  $scratchLoc$'all_cases_dupes'$pheN #  72217
numAllCases=$(wc -l < $scratchLoc$'all_cases'$pheN)

# duplicate this, in case we will need to refer back to the original number of cases (including the ones we may not select)
cp $scratchLoc$'all_cases'$pheN $scratchLoc$'all_cases'$pheN$'_toomany'

# if this is more than 50% than the total (of 125K), then we remove enough to make just 50% )
if [ "$numAllCases" -gt "$maxAllowedCases" ]; then
echo "too many cases"

#  prioritise those that have less missing covars
awk 'FNR == NR { phenos[$2] = $2; next; }
FNR <= NR {  if( $2 in phenos ) { print $0 } }
' $scratchLoc$'all_cases'$pheN$'_toomany' $rawLoc$'indisWithMaxMissing_25' > $rawLoc$'indisWithMaxMissing_25_cases'$pheN
head $rawLoc$'indisWithMaxMissing_25_cases'$pheN
wc -l $rawLoc$'indisWithMaxMissing_25_cases'$pheN # 79150


# intersect from the QC passed ones 
awk 'FNR == NR { phenos[$2] = $2; next; }
FNR <= NR {  if( $2 in phenos ) { print $1"\t"$2 } }
' $rawLoc$'indisWithMaxMissing_25_cases'$pheN $UKBB_PLINK1$'ALL_KEEP.fam' > $rawLoc$'indisWithMaxMissing_25_cases'$pheN$'_QC'
head $rawLoc$'indisWithMaxMissing_25_cases'$pheN$'_QC'
wc -l $rawLoc$'indisWithMaxMissing_25_cases'$pheN$'_QC' #


numAllCases2=$(wc -l < $rawLoc$'indisWithMaxMissing_25_cases'$pheN$'_QC')
echo $numAllCases2

# if this is STILL too many 
if [ "$numAllCases2" -gt "$maxAllowedCases" ]; then
echo "still too many cases"
# we randomly pick from these the correct number via reproducible random number gen
arguments='/home/mk907/scripts/R/randomNums.R '$numAllCases2$' 42 '$scratchLoc$'randomNums_cases'
Rscript $arguments

awk -v maxAllowedCases="$maxAllowedCases" '{if (FNR <= maxAllowedCases) print $0}' $scratchLoc$'randomNums_cases' > $scratchLoc$'randomNums_cases_selected'
head $scratchLoc$'randomNums_cases_selected'

# pick the actual sample ids 
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $0 } }
' $scratchLoc$'randomNums_cases_selected' $rawLoc$'indisWithMaxMissing_25_cases'$pheN$'_QC' > $scratchLoc$'all_cases'$pheN
head $scratchLoc$'all_cases'$pheN

else # otherwise just use however many we found
echo "OK number of cases after QC"
# overwrite the cases 
cp $rawLoc$'indisWithMaxMissing_25_cases'$pheN$'_QC' $scratchLoc$'all_cases'$pheN
fi


else # initial else for already right number of cases
echo "OK number of cases initallly"
fi # end of if there were too many cases first initially



######################
# however we got the correct number of cases list,  we pick now the controls too:
numAllCases=$(wc -l < $scratchLoc$'all_cases'$pheN)
# we make up this with enough controls, so that we would have 100K in the train/valid set, ie 
# 0.8 * X = 100 -> X = 125, 
# IE will need 125K indis in total
numControlsToPick=$(awk -v numAllCases="$numAllCases" -v maxSamplesize="$maxSamplesize" 'BEGIN {print  maxSamplesize- numAllCases}')
echo $numControlsToPick  # 65249


# remove all the cases from the 'low missing' indis (even the ones we did NOT select, as they had too many missing etc!)
awk 'FNR == NR { phenos[$2] = $2; next; }
FNR <= NR {  if( $2 in phenos == 0) { print $0 } }
' $scratchLoc$'all_cases'$pheN$'_toomany' $rawLoc$'indisWithMaxMissing_25' > $rawLoc$'indisWithMaxMissing_25_controls'$pheN
head $rawLoc$'indisWithMaxMissing_25_controls'$pheN
wc -l $rawLoc$'indisWithMaxMissing_25_controls'$pheN # 208328


# intersect from the QC passed ones 
awk 'FNR == NR { phenos[$2] = $2; next; }
FNR <= NR {  if( $2 in phenos ) { print $1"\t"$2 } }
' $rawLoc$'indisWithMaxMissing_25_controls'$pheN $UKBB_PLINK1$'ALL_KEEP.fam' > $rawLoc$'indisWithMaxMissing_25_controls'$pheN'_QC'
head $rawLoc$'indisWithMaxMissing_25_controls'$pheN'_QC'
wc -l $rawLoc$'indisWithMaxMissing_25_controls'$pheN'_QC' # 142671

# from this pick a random $numControlsToPick number of indis 
numTotalIndis=$(wc -l < $rawLoc$'indisWithMaxMissing_25_controls'$pheN'_QC')
# reproducible random number gen
arguments='/home/mk907/scripts/R/randomNums.R '$numTotalIndis$' 42 '$scratchLoc$'randomNums_controls'
Rscript $arguments

awk -v numControlsToPick="$numControlsToPick" '{if (FNR <= numControlsToPick) print $0}' $scratchLoc$'randomNums_controls' > $scratchLoc$'randomNums_controls_selected'
head $scratchLoc$'randomNums_controls_selected'


# pick the actual sample ids 
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $0 } }
' $scratchLoc$'randomNums_controls_selected' $rawLoc$'indisWithMaxMissing_25_controls'$pheN'_QC' > $scratchLoc$'all_controls'$pheN
head $scratchLoc$'all_controls'$pheN


# merge cases and controls create final subset of indis
cat $scratchLoc$'all_cases'$pheN $scratchLoc$'all_controls'$pheN > $rawLoc$'REDUCED_KEEPLIST'$pheN
head $rawLoc$'REDUCED_KEEPLIST'$pheN
wc -l $rawLoc$'REDUCED_KEEPLIST'$pheN # 125000



# pick the actual sample ids 
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $0 } }
' $scratchLoc$'randomNums_training' $rawLoc$'REDUCED_KEEPLIST'$pheN > $rawLoc$'training_keep'$pheN


awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $0 } }
' $scratchLoc$'randomNums_valid' $rawLoc$'REDUCED_KEEPLIST'$pheN > $rawLoc$'valid_keep'$pheN


awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $0 } }
' $scratchLoc$'randomNums_test' $rawLoc$'REDUCED_KEEPLIST'$pheN > $rawLoc$'test_keep'$pheN
head $rawLoc$'test_keep'$pheN
wc -l $rawLoc$'test_keep'$pheN

cat $rawLoc$'training_keep'$pheN $rawLoc$'valid_keep'$pheN > $rawLoc$'trainvalid_keep'$pheN

######
# create PLINK files

# do some basic SNP QC too: the SNPs must be OK in all of the subsets (training/valid/test)
arguments=' --bfile '$UKBB_PLINK1$'ALL_KEEP --memory  43000 --threads 32 --keep '$rawLoc$'training_keep'$pheN$' --exclude '$scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC --maf 0.001 --geno 0.05 --keep-allele-order --write-snplist --out '$scratchLoc$'training_qc'$pheN
/home/mk907/software/plink2/plink2 $arguments

arguments=' --bfile '$UKBB_PLINK1$'ALL_KEEP --memory  43000 --threads 32 --keep '$rawLoc$'valid_keep'$pheN$' --exclude '$scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC --maf 0.001 --geno 0.05 --keep-allele-order --write-snplist --out '$scratchLoc$'valid_qc'$pheN
/home/mk907/software/plink2/plink2 $arguments

arguments=' --bfile '$UKBB_PLINK1$'ALL_KEEP --memory  43000 --threads 32 --keep '$rawLoc$'test_keep'$pheN$' --exclude '$scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC --maf 0.001 --geno 0.05 --keep-allele-order --write-snplist --out '$scratchLoc$'test_qc'$pheN
/home/mk907/software/plink2/plink2 $arguments

# intersect the 3 lists, so that only SNPs are kept that meet the QC criteria in all of the cohorts
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $scratchLoc$'training_qc'$pheN$'.snplist' $scratchLoc$'valid_qc'$pheN$'.snplist' > $scratchLoc$'trainvalid_qc_common'$pheN$'.snplist'
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $scratchLoc$'trainvalid_qc_common'$pheN$'.snplist' $scratchLoc$'test_qc'$pheN$'.snplist' > $scratchLoc$'trainvalidtest_qc_common'$pheN$'.snplist'
head $scratchLoc$'trainvalidtest_qc_common'$pheN$'.snplist'
wc -l $scratchLoc$'trainvalidtest_qc_common'$pheN$'.snplist' # total SNPs


# actually, also create subset for the training+valid too for even smaller RAM footprint
arguments=' --bfile '$UKBB_PLINK1$'ALL_KEEP --memory  43000 --threads 32 --keep '$rawLoc$'trainvalid_keep'$pheN$' --extract '$scratchLoc$'trainvalidtest_qc_common'$pheN$'.snplist --keep-allele-order --make-bed --out '$UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN
/home/mk907/software/plink2/plink2 $arguments

arguments=' --bfile '$UKBB_PLINK1$'ALL_KEEP --memory  43000 --threads 32 --keep '$rawLoc$'test_keep'$pheN$' --extract '$scratchLoc$'trainvalidtest_qc_common'$pheN$'.snplist --keep-allele-order --make-bed --out '$UKBB_PLINK1$'ALL_KEEP_TEST'$pheN
/home/mk907/software/plink2/plink2 $arguments
}

##############################################


# converts my sumstats format into the PRS-CS format, postfixing input with '_PRSCS' (specify if binary)
function convertToPRSCS { 
sumstatsLoc=$1
isBinary=$2

# convert my _sumstats format:
# chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
# to PRS-CS:
# SNP          A1   A2   BETA/OR      P

awk -v isBinary="$isBinary" ' 
FNR <= NR {
if (FNR == 1) {
coefName="BETA"
if(isBinary == "1") {coefName="OR"}
print "SNP\tA1\tA2\t"coefName"\tP" 
}
else { 
COEF=$7
if(isBinary == "1") {COEF=exp($7)}
print $3"\t"$4"\t"$5"\t"COEF"\t"$9 
}
}
' OFS='\t' $sumstatsLoc > $sumstatsLoc$'_PRSCS'

#head $sumstatsLoc$'_PRSCS'

}
export -f convertToPRSCS # this makes local functions executable when bsubbed


# Performs PRS-CS for 1 population (EUR)
eursums=$rawPRSLoc$'_PRSCS'
Neur=$Neur
validBim=$validBim
pheno=$pheno_name
chrom=$i
ncores=16
scratchL=$scratchLoca
function performPRSCS { 
eursums=$1
Neur=$2
validBim=$3
pheno=$4
chrom=$5
ncores=$6
scratchL=$7

export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

echo "Made directory at "$scratchL$pheno$"/"
prscsxdir="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/prscsscript/PRScsx/"
prscsrefs="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/prscsrefs/"
#scratchLocV4="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/scratch/"
mkdir -p $scratchL$pheno$"/"
/usr/local/software/master/miniconda/2/bin/python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$validBim --sst_file=$eursums --n_gwas=$Neur --pop=EUR --out_dir=$scratchL$pheno/ --out_name=$pheno --chrom=$chrom  --seed=42

}
export -f performPRSCS # this makes local functions executable when bsubbed




# rawPRSLoc=$ukbbres$'_sumstats'
# pheno_name=$phename
# scratchLoca=$scratchLoc
# isBina='0'
# performs PRS-CS for each of the 22 chroms
# b=0
# order=1
# rawPRSLoc=$ukbbres$'_sumstats'
# pheno_name=$phename
# scratchLoca=$scratchLoc 
# isBina='0'
# i=22



function PRSCS_allChroms { 
rawPRSLoc=$1
pheno_name=$2
scratchLoca=$3
isBina=$4
convertToPRSCS $rawPRSLoc $isBina


# get AVERAGE sample size, this was better than minimum
Neur=$(awk '{ if(FNR > 1) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawPRSLoc)
echo $Neur


validBim=$rawPRSLoc$'_validBim'
# create a 'fake' .bim file for PRSx from my sumstats file that has all the info
awk '{if (FNR !=1) {print $1"\t"$3"\t0\t"$2"\t"$4"\t"$5} }' $rawPRSLoc > $validBim$'.bim'


for ((i=1; i<=22; i++)); do # $numChroms
mkdir -p $scratchLoca$pheno_name$'/'
mkdir -p $scratchLoca$'logs/'
# 1) PRS-CS
# check if output doesn't already exist
outfile=$scratchLoca$pheno_name$'/'$pheno_name$'_EUR_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
bq=1
else
bq=1
echo "PRS-CS SUBMITTING: "$pheno_name$' '$i

sbatch --mem 43000 --cpus-per-task 16 --time 4:0:0 --partition cclake --account danesh-sl3-cpu -e $scratchLoca$'logs/'$pheno_name$'_'$i$'.err' -o $scratchLoca$'logs/'$pheno_name$'_'$i$'.out' --wrap "performPRSCS $rawPRSLoc$'_PRSCS' $Neur $validBim $pheno_name $i 16 $scratchLoca"
#performPRSCS $rawPRSLoc$'_PRSCS' $Neur $validBim $pheno_name $i 16 $scratchLoca

fi

done # end of chrom loop

} 
export -f PRSCS_allChroms # this makes local functions executable when bsubbed


# Builds the individual profile scores for all chromosomes for PRS-CS
function BuildPRS_allChroms_PRSCS { 
penoNam=$1
scratchLoca=$2
bfile=$3
datLoc=$4
testPheno=$5

#  concat all chrom PRS into single files
cat $scratchLoca$penoNam$'/'$penoNam$'_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $scratchLoca$penoNam$'/'$penoNam$'_PRSCS'
head $scratchLoca$penoNam$'/'$penoNam$'_PRSCS'
wc -l $scratchLoca$penoNam$'/'$penoNam$'_PRSCS' 

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $scratchLoca$penoNam$'/'$penoNam$'_PRSCS' > $scratchLoca$penoNam$'/'$penoNam$'_PRSCS_no_dupes'

BuildPRS_cluster_PRSCS $scratchLoca$penoNam$'/'$penoNam$'_PRSCS_no_dupes' $penoNam$'_PRSCS' $testPheno '1' $bfile $datLoc

}
export -f BuildPRS_allChroms_PRSCS # this makes local functions executable when bsubbed





# PRS_rawLoc=$scratchLoca$penoNam$'/'$penoNam$'_PRSCS_no_dupes' 
# phe_name=$penoNam$'_PRSCS' 
# testPheno=$testPheno 
# PRSCS='1' 
# bfil=$bfile 
# datLo=datLoc
# builds a PRS from PRS-CS or LDPred2 formatted score file
function BuildPRS_cluster_PRSCS {
PRS_rawLoc=$1
phe_name=$2
testPheno=$3
PRSCS=$4
bfil=$5
datLo=$6


i=21

# generate PRS per each chrom
mkdir -p $datLo$'PRSProfiles/'
mkdir -p $datLo$'PRSProfiles/'$phe_name$'/'
for ((i=1; i<=22; i++)); do

# only submit if it doesn't exist yet
outfile=$datLo$'PRSProfiles/'$phe_name$'/'$i$'.sscore'
rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
bqq=1
else
bqq=1

#PRS-CS
if [ $PRSCS == "1" ] ; then
echo "PRS-CS PRS"
arguments=' --bfile '$bfil$' --memory  43000 --threads 16 --score '$PRS_rawLoc$' 2 4 6 ignore-dup-ids cols=-scoreavgs,+scoresums --pheno '$testPheno$' --out '$datLo$'PRSProfiles/'$phe_name$'/'$i

else
# LDpred2/Rapido
echo "LDpred2 PRS"
arguments=' --bfile '$bfil$' --memory  43000 --threads 16 --score '$PRS_rawLoc$' ignore-dup-ids cols=-scoreavgs,+scoresums --pheno '$testPheno$' --out '$datLo$'PRSProfiles/'$phe_name$'/'$i
fi 

#sbatch --mem 43000 --cpus-per-task 16 --time 4:0:0 --partition cclake --account danesh-sl3-cpu --wrap "/home/mk907/software/plink2/plink2 $arguments"
$plink2 $arguments

fi # if score exists

done

}
export -f BuildPRS_cluster_PRSCS # this makes local functions executable when bsubbed









# sums each 22 chroms PRS into a single profile score, is not tied to V4
function SumPRS_Agnostic {
phe_name=$1
datLoc=$2


i=21

# sum across each chrom
for ((i=1; i<=22; i++)); do

# if it is first iteration, we just copy the score file for 1st chrom, as 
if [ $i == 1 ] ; then
echo "first chrom"
# 3rd col is the true pheno, 6th col is the scoresum, we only want a IID,pheno,PRS file
                        
awk '{print $2"\t"$3"\t"$6 }' $datLoc$'PRSProfiles/'$phe_name$'/'$i$'.sscore' > $datLoc$'PRSProfiles/'$phe_name$'_all.sscore'
cp $datLoc$'PRSProfiles/'$phe_name$'_all.sscore' $datLoc$'PRSProfiles/'$phe_name$'_all.temp'

else

awk 'FNR == NR { sums[ $1 ] = $3; next; }
FNR <= NR { if( $2 in sums ) {print $2"\t"$3"\t"$6+sums[$2] } }
' $datLoc$'PRSProfiles/'$phe_name$'_all.temp' $datLoc$'PRSProfiles/'$phe_name$'/'$i$'.sscore' > $datLoc$'PRSProfiles/'$phe_name$'_all.sscore'
cp $datLoc$'PRSProfiles/'$phe_name$'_all.sscore' $datLoc$'PRSProfiles/'$phe_name$'_all.temp'

fi 

done

rm -rf $datLoc$'PRSProfiles/'$phe_name$'_all.temp'


}
export -f SumPRS_Agnostic # this makes local functions executable when bsubbed




# Quantitative traits: converts the results of a new version of PLINK2 to the common sumstats format
function ConvertPLINK2QuantitativePhenoToSumtats_new { 
numIndis=$1 
plink2file=$2
outlocation=$3


#   1 		2     3              4       5       6                      7       8       9               10      11      12      13      14      15      16      17      18
# CHROM    POS    ID             REF     ALT     PROVISIONAL_REF?       A1      OMITTED A1_FREQ         TEST    OBS_CT  BETA    SE      L95     U95     T_STAT  P       ERRCODE
#1       754182  rs3131969       G       A       Y                      A       G       0.128498        ADD     74433   0.00677307      0.00775579      -0.00842801     0.0219741 0.873291        0.382507 

# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
#  1     2     3    4     5        6         7    8   9   10

awk  -v numIndis=$numIndis ' FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { 
# find the non-effect allele, this could either be REF or ALT
other_allele=$4 # default to REF
if(other_allele == $7) {other_allele = $5} # if REF is the same as A1, then use ALT

print $1"\t"$2"\t"$3"\t"$7"\t"other_allele"\t"$9"\t"$12"\t"$13"\t"$17"\t"numIndis
}
} ' OFS='\t' $plink2file > $outlocation
}
export -f ConvertPLINK2QuantitativePhenoToSumtats_new # this makes local functions executable when bsubbed




