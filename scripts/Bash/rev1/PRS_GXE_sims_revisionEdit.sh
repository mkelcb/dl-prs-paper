

# Revision 1 for Sims bash script
# assumes PRS_GX_functions_v2 loaded in, and PRS_GXE_sims.sh already performed.

# repeats analysis with adding an additive+epistatic scenario

# screen -r -D 3115212.pts-202.login-q-1


#####################

source ~/pytorch-env113/bin/activate

########
# take the _smallcases genotype set, randomly subset to 500K SNPs
pheN="_smallcases"


$simsData


#############################################
# Create mixed additive+epistatic scenarios, by combining the additive and 4-way interaction phenos
#############################################


# mixes two phenotpes, the additive and epistatic scenarios, with given proportions
function mixPhenos { 
additive=$1
epistatic=$2
additive_ratio=$3
epistatic_ratio=$4
out=$5

awk -v additive_ratio="$additive_ratio" -v epistatic_ratio="$epistatic_ratio"  'FNR == NR { file1[ $1 ] = $3;next; } 
FNR <= NR {  if( $1 in file1) {print $1"\t"$2"\t"additive_ratio*file1[$1]+epistatic_ratio *$3 } }
'  FS="\t" $additive  $epistatic > $out

}






order=4
b=19
for ((b=0; b<numSimReplicates; b++)); do 

# create 2 mixes between additive and epistatic scenarios with the following ratios: 1:2 and 2:1


# 1:2
mixPhenos $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_1_v'$b$'.pheno' $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_4_v'$b$'.pheno' 1 2 $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_mixed12_v'$b$'.pheno'
mixPhenos $simsDataTruthset$'phenos/GWAS_KNET_test_1_v'$b$'.pheno' $simsDataTruthset$'phenos/GWAS_KNET_test_4_v'$b$'.pheno' 1 2 $simsDataTruthset$'phenos/GWAS_KNET_test_mixed12_v'$b$'.pheno'

# 2:1
mixPhenos $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_1_v'$b$'.pheno' $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_4_v'$b$'.pheno' 2 1 $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_mixed21_v'$b$'.pheno'
mixPhenos $simsDataTruthset$'phenos/GWAS_KNET_test_1_v'$b$'.pheno' $simsDataTruthset$'phenos/GWAS_KNET_test_4_v'$b$'.pheno' 2 1 $simsDataTruthset$'phenos/GWAS_KNET_test_mixed21_v'$b$'.pheno'


# $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno'
# $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno'

done # end of numVersions  loop



$simsDataTruthset

# need to create pheno file that has IIDs and pheno in 2 cols
$simsDataTruthset$"phenos/"


# 1:2
$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_mixed12_v'$b$'.pheno'
$simsDataTruthset$'phenos/GWAS_KNET_test_mixed12_v'$b$'.pheno'

# 2:1
$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_mixed21_v'$b$'.pheno'
$simsDataTruthset$'phenos/GWAS_KNET_test_mixed21_v'$b$'.pheno'



#############################################
# III. perform basic PLINK GWAS + PRS-CS
#############################################
N_eff_UKBB=$(wc -l < $rawLoc$'training_keep_smallcases')
b=1
order=1



mixedScenarios=( 'mixed12' 'mixed21' )
arraylength=${#mixedScenarios[@]}
for ((b=0; b<numSimReplicates; b++)); do 
for (( j=1; j<${arraylength}+1; j++ )); do
order=${mixedScenarios[$j-1]}

echo "order is: "$order

phename='GWAS_'$order$'_v'$b
ukbbres=$resultsLoc$'sims/'$phename$'.PHENO1.glm.linear'

# cp $ukbbres $resultsLoc$'sims/'$ukbbres$'_old'
# GWAS is restricted to the training set, as the PRS weights will be fed back to the DL model via the pre-multiplication thing...
arguments=' --bfile '$UKBB_PLINK1$'SIMS_TRAINVALID --memory  40000 --threads 16  --linear hide-covar allow-no-covars cols=+err,+a1freq  --mac 20 --keep '$rawLoc$'training_keep_smallcases --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno --out '$resultsLoc$'sims/'$phename$' --ci 0.95 --allow-extra-chr --snps-only --bp-space 1'
#sbatch --mem 40000 --cpus-per-task 16 --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$plink2 $arguments"
$plink2 $arguments

# convert to my sumstats format
ConvertPLINK2QuantitativePhenoToSumtats_new $N_eff_UKBB $ukbbres $ukbbres$'_sumstats'


# perform PRS-CS
PRSCS_allChroms $ukbbres$'_sumstats' $phename $scratchLoc '0'



done  # end of numOrders  loop

done # end of numVersions  loop



#########
# Evaluate PRS-CS PRS
# test set phenos

mixedScenarios=( 'mixed12' 'mixed21' )
arraylength=${#mixedScenarios[@]}

b=0
for ((b=0; b<numSimReplicates; b++)); do 
for (( j=1; j<${arraylength}+1; j++ )); do
order=${mixedScenarios[$j-1]}

echo "order is: "$order
phename='GWAS_'$order$'_v'$b

# Build PRS
BuildPRS_allChroms_PRSCS $phename $scratchLoc $UKBB_PLINK1$'SIMS_TEST' $simsData $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno'

# Sum and Evaluate
SumPRS_Agnostic $phename$'_PRSCS' $simsData

Evaluate_PRS $simsData$'PRSProfiles/'$phename$'_PRSCS_all.sscore' $resultsLoc$'sims/GWAS_'$order$'_res' '0' '1' # append across ALL replicates!

# also need to convert the $scratchLoc$phename$'/'$phename$'_PRSCS_no_dupes' into my Knet per SNP PRS format:
# map the headerless:
# 10      rs11591988      126070  T       C       1.161153e-05
# to 
#hm_chr  hm_pos  effect_allele   other_allele    effect_weight
#1       754182  G               A               0.0001081757
awk '{if (FNR == 1) {print "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight"} print $1"\t"$3"\t"$4"\t"$5"\t"$6}' $scratchLoc$phename$'/'$phename$'_PRSCS_no_dupes' > $simsData$phename$knetPRSext


done  # end of numOrders  loop

done # end of numVersions  loop







############################################## 
# V. Sim Analyses (to be run on the GPU nodes) 
############################################## 


mixedScenarios=( 'mixed12' 'mixed21' )
arraylength=${#mixedScenarios[@]}
b=19
# fit KNET, in 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for (( j=1; j<${arraylength}+1; j++ )); do
order=${mixedScenarios[$j-1]}

phename='GWAS_'$order$'_v'$b
echo "order is: "$order

# mixed scenarios
mkdir -p $resultsLoc$'sims/'$phename

outName='res_'$order
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments  > $resultsLoc$'sims/'$phename$'/'$outName$'_args'


done  # end of numOrders  loop

done # end of numVersions  loop




mixedScenarios=( 'mixed12' 'mixed21' )
arraylength=${#mixedScenarios[@]}
# check if all results are there...
for ((b=0; b<numSimReplicates; b++)); do 
for (( j=1; j<${arraylength}+1; j++ )); do
order=${mixedScenarios[$j-1]}
phename='GWAS_'$order$'_v'$b

# mixed scenarios
outName='res_'$order
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'

done  # end of numOrders  loop

done # end of numVersions  loop




####################
# Inference KNET:
mixedScenarios=( 'mixed12' 'mixed21' )
arraylength=${#mixedScenarios[@]}
b=0
phename='GWAS_'$order$'_v'$b
#  KNET, in the same 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for (( j=1; j<${arraylength}+1; j++ )); do
order=${mixedScenarios[$j-1]}
phename='GWAS_'$order$'_v'$b

# mixed scenarios
outName='res_'$order
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"


done  # end of numOrders  loop
done # end of numVersions  loop



# Check if all results are there...
mixedScenarios=( 'mixed12' 'mixed21' )
arraylength=${#mixedScenarios[@]}
#  KNET, in the same 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for (( j=1; j<${arraylength}+1; j++ )); do
order=${mixedScenarios[$j-1]}
phename='GWAS_'$order$'_v'$b

echo "order is: "$order
# mixed scenarios
outName='res_'$order
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName

done  # end of numOrders  loop
done # end of numVersions  loop

###############





mixedScenarios=( 'mixed12' 'mixed21' )
arraylength=${#mixedScenarios[@]}
phename='GWAS_'$order$'_v'$b
mkdir -p $resultsLoc$phename$'/eval/'
#  KNET, in the same 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for (( j=1; j<${arraylength}+1; j++ )); do
order=${mixedScenarios[$j-1]}
phename='GWAS_'$order$'_v'$b
mkdir -p $resultsLoc$phename$'/eval/'


echo "order is: "$order
# mixed scenario
outName='res_'$order
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'


done  # end of numOrders  loop
done # end of numVersions  loop




#############


# extract all results into 3 text files ( we work off the correlation not correlation^2)
header="replicate\tres_12\tres_21\tPRSIndi_12\tPRSIndi_21"
echo -e $header > $resultsLoc$'simResRev_nonlinear'
echo -e $header > $resultsLoc$'simResRev_noAct'
echo -e $header > $resultsLoc$'simResRev_linear'

for ((b=0; b<numSimReplicates; b++)); do 



# 1. Mixed 1:2
order='mixed12'
phename='GWAS_'$order$'_v'$b
outName='res_'$order
res_12_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_12_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_12_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# 2. Mixed 2:1
order='mixed21'
phename='GWAS_'$order$'_v'$b
outName='res_'$order
res_21_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_21_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_21_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# get the linear PRS:
PRSIndi_12=$(awk -v b="$b" '{ if(FNR == (b+1) ) {print $2 } }' $resultsLoc$'sims/GWAS_mixed12_res_r2_redux')
PRSIndi_21=$(awk -v b="$b" '{ if(FNR == (b+1) ) {print $2 } }' $resultsLoc$'sims/GWAS_mixed21_res_r2_redux')




echo -e $b$"\t"$res_12_nonlinear$"\t"$res_21_nonlinear$"\t"$PRSIndi_12$"\t"$PRSIndi_21 >> $resultsLoc$'simResRev_nonlinear'
echo -e $b$"\t"$res_12_linear$"\t"$res_21_linear$"\t"$PRSIndi_12$"\t"$PRSIndi_21 >> $resultsLoc$'simResRev_linear'
echo -e $b$"\t"$res_12_noAct$"\t"$res_21_noAct$"\t"$PRSIndi_12$"\t"$PRSIndi_21 >> $resultsLoc$'simResRev_noAct'

done # end of numVersions  loop

#!!! HERE

#############################################
# VI. Plot Results
#############################################


#$resultsLoc$'simResRev_nonlinear'
#$resultsLoc$'simResRev_noAct'
#$resultsLoc$'simResRev_linear'

arguments='/home/mk907/scripts/R/PRS_GXE_SimsPlotterRev.R '$resultsLoc$'simResRev '$resultsLoc$'simResRev_res PRS_GxE_simulationsRev'
Rscript $arguments
head $resultsLoc$'simResRev_res_tab'



##results: 
# scenario	            mean_nonlinear	        t-test_p
# res_21                0.0281966389116897      0.000110870465792107
# res_12	            0.0937996699002529      5.71999191257027e-05

# a significant positive value means that the non-linear model is better

# Expectation:

# res_21: additive:epistatic: 66-33%: nonlinear should be worst, as only 1/3 is epistatic
# res_12: additive:epistatic: 33-66%: nonlinear should be closest to 100% epistatic, as most of the phenotype should come from nonlinrear effects

# res_12: good as it is 9.3%, close to the 9.5% when epistasis is 100%
# res_21: OK, as it is 3%, which is appropriately lower and in line with the expectations  



