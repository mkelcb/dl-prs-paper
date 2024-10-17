

# PRS GXE project simulations
# 

# 

#####################

source ~/pytorch-env113/bin/activate

########
# take the _smallcases genotype set, randomly subset to 500K SNPs
pheN="_smallcases"


mkdir -p $simsData
# from this pick a random $numControlsToPick number of indis 
numSNPs=$(wc -l < $UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$'.bim')
# reproducible random number gen
arguments='/home/mk907/scripts/R/randomNums.R '$numSNPs$' 42 '$scratchLoc$'randomNums_SNPs'
Rscript $arguments

awk -v numSNPsToPick="$numSNPsToPick" '{if (FNR <= numSNPsToPick) print $0}' $scratchLoc$'randomNums_SNPs' > $scratchLoc$'randomNums_SNPs_selected'
head $scratchLoc$'randomNums_SNPs_selected'


# pick the actual  SNPs
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $2 } }
' $scratchLoc$'randomNums_SNPs_selected' $UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$'.bim' > $scratchLoc$'all_SNPs'
head $scratchLoc$'all_SNPs'




arguments=' --bfile '$UKBB_PLINK1$'ALL_KEEP_TRAINVALID'$pheN$' --memory  43000 --threads 32 --extract '$scratchLoc$'all_SNPs --keep-allele-order --make-bed --out '$UKBB_PLINK1$'SIMS_TRAINVALID'
/home/mk907/software/plink2/plink2 $arguments

arguments=' --bfile '$UKBB_PLINK1$'ALL_KEEP_TEST'$pheN$' --memory  43000 --threads 32 --extract '$scratchLoc$'all_SNPs --keep-allele-order --make-bed --out '$UKBB_PLINK1$'SIMS_TEST'
/home/mk907/software/plink2/plink2 $arguments



#############################################
# II. Simulate 20 phenos with 0.5 h2, from 2000 causal vars, from 4 way interactions or additive
#############################################
# 500000 * 0.004 =  2000

mkdir -p $simsDataTruthset

arguments=' --out '$simsDataTruthset$'sim --plinkFile1 '$UKBB_PLINK1$'SIMS_TRAINVALID --plinkFile2 '$UKBB_PLINK1$'SIMS_TEST --percCausal '$causalPerc$' --h2 '$h2level$' --numVersions '$numSimReplicates$' --numInteractions -1'
#python3 /home/mk907/software/knet2/phenoSim2.py $arguments
sbatch --mem 128000 --cpus-per-task 3 --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "python3 /home/mk907/software/knet2/phenoSim2.py $arguments"






# need to create pheno file that has IIDs and pheno in 2 cols
mkdir -p $simsDataTruthset$"phenos/"
awk '{ {print $2"\t"$2}}' $UKBB_PLINK1$'SIMS_TRAINVALID.fam' > $scratchLoc"GWAS_KNET_100_trainvalid_indis"
awk '{ {print $2"\t"$2}}' $UKBB_PLINK1$'SIMS_TEST.fam' > $scratchLoc"GWAS_KNET_100_test_indis"

order=4
b=19
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 


paste $scratchLoc"GWAS_KNET_100_trainvalid_indis" $simsDataTruthset$'sim_train_'$order$'_v'$b > $simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno'
paste $scratchLoc"GWAS_KNET_100_test_indis" $simsDataTruthset$'sim_valid_'$order$'_v'$b > $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno'



done  # end of numOrders  loop

done # end of numVersions  loop

#############################################
# III. perform basic PLINK GWAS + PRS-CS
#############################################
N_eff_UKBB=$(wc -l < $rawLoc$'training_keep_smallcases')
b=1
order=1
mkdir -p $resultsLoc$'sims/'
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 

# only do linear and 4th order interaction ones
if [ "$order" = "1" ] || [ "$order" = "4" ]; then
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

fi


done  # end of numOrders  loop

done # end of numVersions  loop





#########
# Evaluate PRS-CS PRS
# test set phenos
$UKBB_PLINK1$'SIMS_TEST'
$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno'
rm -rf $resultsLoc$'sims/GWAS_1_res_r2'
rm -rf $resultsLoc$'sims/GWAS_4_res_r2'
rm -rf $resultsLoc$'sims/GWAS_1_res_r2_redux'
rm -rf $resultsLoc$'sims/GWAS_4_res_r2_redux'

b=0
order=1
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 

# only do linear and 4th order interaction ones
if [ "$order" = "1" ] || [ "$order" = "4" ]; then
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
fi

done  # end of numOrders  loop

done # end of numVersions  loop



#############################################
# IV. generate Haplotype effect scenarios
#############################################

# create r2 of all causal variants against all others
# for each additive simulation scenario, and calculate r2 of the causal SNPs against all else
# find all which had 2+ tagging variants with r2 >= 0.25
order=1
b=0
for ((b=0; b<numSimReplicates; b++)); do 
phename='GWAS_'$order$'_v'$b
ukbbres=$resultsLoc$'sims/GWAS_'$order$'_v'$b$'.PHENO1.glm.linear'

# find the causal SNPs, these were saved as indices by the simulation script, which refer to rownumbers in the .bim file, so match these back
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $2 } }
' $simsDataTruthset$'sim_1_v'$b$'.indices' $UKBB_PLINK1$'SIMS_TRAINVALID.bim' > $scratchLoc$'sim_1_v'$b$'.SNP'
head $scratchLoc$'sim_1_v'$b$'.SNP'

arguments=' --memory 40000 --bfile '$UKBB_PLINK1$'SIMS_TRAINVALID  --tag-r2 0.5 --list-all --tag-kb 500 --show-tags '$scratchLoc$'sim_1_v'$b$'.SNP --out '$scratchLoc$'sim_1_v'$b$'_tags --keep-allele-order --allow-extra-chr --allow-no-sex '
$plink $arguments



# find the SNPs that have at least 2 tagging SNPs 
# $scratchLoc$'sim_1_v'$b$'_tags.tags.list' has signature:
#      SNP  CHR         BP NTAG       LEFT      RIGHT   KBSPAN TAGS
# rs2542334   22   16694612    2   16693517   16695440    1.923 rs415170|rs2587108
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { 
# if it has NTAG 2 or more
if(FNR > 1 && $4 >= 2) {
# split the rs415170|rs2587108 into array of rs415170,rs2587108...
split($8,a,"|"); 
numtaggingSNPs=0
# go through each SNP, but only count those that are NOT themselves targets
for (snp in a) {
if ( snp in test == 0 ) {numtaggingSNPs++}
}
if(numtaggingSNPs >= 2) {print $1}
} }
' $scratchLoc$'sim_1_v'$b$'.SNP' $scratchLoc$'sim_1_v'$b$'_tags.tags.list' > $scratchLoc$'sim_1_v'$b$'.SNP.excludable'
head $scratchLoc$'sim_1_v'$b$'.SNP.excludable'
wc -l  $scratchLoc$'sim_1_v'$b$'.SNP.excludable'

# exclude a random 50% of these excludable SNPs 
# from this pick a random $numControlsToPick number of indis 
numSNPs=$(wc -l < $scratchLoc$'sim_1_v'$b$'.SNP.excludable')

# if we have more than the half (1000) we exclude only up to the half, but if we have less, than we exclude whatever we can
num_half=$(awk -v numSNPs=$numSNPs -v halfCausals=$halfCausals 'BEGIN{ if (halfCausals > numSNPs) {print numSNPs} else {print halfCausals}  }')
echo "excluding half ("$num_half$") of total numSNPs: "$numSNPs



# reproducible random number gen
arguments='/home/mk907/scripts/R/randomNums.R '$numSNPs$' 42 '$scratchLoc$'randomNums_SNPs'
Rscript $arguments

awk -v num_half="$num_half" '{if (FNR <= num_half) print $0}' $scratchLoc$'randomNums_SNPs' > $scratchLoc$'randomNums_SNPs_selected'
head $scratchLoc$'randomNums_SNPs_selected'

# pick the actual SNPs
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $2 } }
' $scratchLoc$'randomNums_SNPs_selected' $UKBB_PLINK1$'SIMS_TRAINVALID.bim' > $scratchLoc$phename'excluded_SNPs'
head $scratchLoc$phename'excluded_SNPs'

# create PLINK files (must do this as KNET will only work on PLINK files, cannot work on subsets...)
arguments=' --bfile '$UKBB_PLINK1$'SIMS_TRAINVALID --memory  43000 --threads 32 --exclude '$scratchLoc$phename'excluded_SNPs --keep-allele-order --make-bed --out '$UKBB_PLINK1$'SIMS_TRAINVALID_haplo_'$b
/home/mk907/software/plink2/plink2 $arguments

arguments=' --bfile '$UKBB_PLINK1$'SIMS_TEST --memory  43000 --threads 32 --exclude '$scratchLoc$phename'excluded_SNPs --keep-allele-order --make-bed --out '$UKBB_PLINK1$'SIMS_TEST_haplo_'$b
/home/mk907/software/plink2/plink2 $arguments



# Filter out 
# 1. create clumps (default)

rm -rf $resultsLoc$'sims/'$phename$'_Clump.clumps'
arguments=' --bfile '$UKBB_PLINK1$'SIMS_TRAINVALID --memory  43000 --threads 32 --clump '$ukbbres$' --clump-kb 500 --clump-r2 0.10  --out '$resultsLoc$'sims/'$phename$'_Clump  --allow-extra-chr --snps-only --bp-space 1 --allow-no-sex'
$plink2 $arguments
rm -rf $resultsLoc$'sims/'$phename$'_Clump.log'

#date -r $UKBB_PLINK1$'SIMS_TEST_haplo_'$b$'.fam'
# $resultsLoc$'sims/'$phename$'_Clump.clumps'
# 2. thin clumps so that no 2 clumps are closer than 500kb (always keep the one with lower p-val)
# (we need this as ld pruning /clumping checks for ld between the EXISTING panel of SNPs, whereas we are worried about LD with missing SNPs NOT on the panel, whereas the 2 SNPs in the panel could be in perfect LE, which could impute missing causal variants via interactions)
arguments='/home/mk907/scripts/R/clumpThinner.R '$resultsLoc$'sims/'$phename$'_Clump.clumps '$resultsLoc$'sims/'$phename$'_Clump_kept 500'
Rscript $arguments

arguments=' --memory 43000 --bfile '$UKBB_PLINK1$'SIMS_TRAINVALID --threads 16 --allow-no-sex --extract '$resultsLoc$'sims/'$phename$'_Clump_kept --keep-allele-order --make-bed --out '$UKBB_PLINK1$'SIMS_TRAINVALID_haplo_filtered'$b
/home/mk907/software/plink2/plink2 $arguments

arguments=' --memory 43000 --bfile '$UKBB_PLINK1$'SIMS_TEST --threads 16 --allow-no-sex --extract '$resultsLoc$'sims/'$phename$'_Clump_kept --keep-allele-order --make-bed --out '$UKBB_PLINK1$'SIMS_TEST_haplo_filtered'$b
/home/mk907/software/plink2/plink2 $arguments

done # end of numVersions  loop



#####################


# create PRS-CS for the haplotype effect ones too (but NOT for the filtered ones, as with those we eliminated the haplo effect via filtering)
order=1 # only do it for additives
b=0
for ((b=0; b<numSimReplicates; b++)); do 
phename='GWAS_'$order$'_v'$b
ukbbres=$resultsLoc$'sims/'$phename$'.PHENO1.glm.linear'

# remove the variants from the summary statistics
awk 'FNR == NR { test[ $2 ] = $2; next; } FNR <= NR { if( FNR ==1 || $3 in test ) {print $0 } }
' $UKBB_PLINK1$'SIMS_TRAINVALID_haplo_'$b$'.bim' $ukbbres$'_sumstats' > $ukbbres$'_sumstats_haplo'

# perform PRS-CS
PRSCS_allChroms $ukbbres$'_sumstats_haplo' $phename$'_haplo' $scratchLoc '0'

done # end of numVersions  loop


#########
# Evaluate PRS-CS PRS
# test set phenos
rm -rf $resultsLoc$'sims/GWAS_1_res_haplo_r2'
rm -rf $resultsLoc$'sims/GWAS_1_res_haplo_r2_redux'

order=1 # only do it for additives
for ((b=0; b<numSimReplicates; b++)); do 
phename='GWAS_'$order$'_v'$b

# Build PRS
BuildPRS_allChroms_PRSCS $phename$'_haplo' $scratchLoc $UKBB_PLINK1$'SIMS_TEST' $simsData $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno'

# Sum and Evaluate
SumPRS_Agnostic $phename$'_haplo_PRSCS' $simsData

Evaluate_PRS $simsData$'PRSProfiles/'$phename$'_haplo_PRSCS_all.sscore' $resultsLoc$'sims/GWAS_'$order$'_res_haplo' '0' '1' # correlation_sq 0.00437 (sd: 1.56e-05 )  AUC 0.69 CI low: 0.6729  / high:  0.7071


# also need to convert the $scratchLoc$phename$'/'$phename$'_PRSCS_no_dupes' into my Knet per SNP PRS format:
# map the headerless:
# 10      rs11591988      126070  T       C       1.161153e-05
# to 
#hm_chr  hm_pos  effect_allele   other_allele    effect_weight
#1       754182  G               A               0.0001081757
awk '{if (FNR == 1) {print "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight"} print $1"\t"$3"\t"$4"\t"$5"\t"$6}' $scratchLoc$phename$'_haplo/'$phename$'_haplo_PRSCS_no_dupes' > $simsData$phename$knetPRSext$'_haplo'

done # end of numVersions  loop





############################################## 
# V. Sim Analyses (to be run on the GPU nodes) 
# (THE BELOW MAY NEED TO BE SUBMITTED IN 2 PARTS, AS THE CONSOLE GETS CHOKED ON THE AMOUNT OF TEXT WHICH THEN RESULTS IN WRONG PARAMS AND ERRORS...)
############################################## 

b=19
order=4
# fit KNET, in 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 
phename='GWAS_'$order$'_v'$b
echo "order is: "$order

# additive scenarios
if [ "$order" = "1" ] ; then
mkdir -p $resultsLoc$'sims/'$phename
# 1. baseline - additive truth (full genotype)
outName='res'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments > $resultsLoc$'sims/'$phename$'/'$outName$'_args'

# 1b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_PRS'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno  --prs '$simsData$phename$knetPRSext$' --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments > $resultsLoc$'sims/'$phename$'/'$outName$'_args'

# Haplotype effect tests (only order 1, ie additive)
# 3a. haplotype-effect - additive truth (incomplete coverage genotype)
outName='res_haplo'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID_haplo_'$b$' --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments > $resultsLoc$'sims/'$phename$'/'$outName$'_args'

# 3b. same as a. but with PRS-CS premultiplied SNPs (incomplete coverage genotype)
outName='res_haplo_PRS'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID_haplo_'$b$' --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno  --prs '$simsData$phename$knetPRSext$'_haplo --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments > $resultsLoc$'sims/'$phename$'/'$outName$'_args'

# 3c. same as a. but with filtered haplotype effects (LD filtered genotype)
outName='res_haplo_filtered'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID_haplo_filtered'$b$' --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments > $resultsLoc$'sims/'$phename$'/'$outName$'_args'
fi

# 4th order interaction scenarios
if [ "$order" = "4" ] ; then
mkdir -p $resultsLoc$'sims/'$phename

# 2. baseline - epistatic truth (full genotype)
outName='res_epi'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments > $resultsLoc$'sims/'$phename$'/'$outName$'_args'

# 2b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_epi_PRS'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TRAINVALID --pheno '$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno  --prs '$simsData$phename$knetPRSext$' --validSet '$rawLoc$'valid_keep_smallcases --saveWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
echo $arguments > $resultsLoc$'sims/'$phename$'/'$outName$'_args'
fi

done  # end of numOrders  loop

done # end of numVersions  loop





# check if all results are there...
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 
phename='GWAS_'$order$'_v'$b

# additive scenarios
if [ "$order" = "1" ] ; then

# 1. baseline - additive truth (full genotype)
outName='res'
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'

# 1b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_PRS'
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'

# Haplotype effect tests (only order 1, ie additive)
# 3a. haplotype-effect - additive truth (incomplete coverage genotype)
outName='res_haplo'
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'

# 3b. same as a. but with PRS-CS premultiplied SNPs (incomplete coverage genotype)
outName='res_haplo_PRS'
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'

# 3c. same as a. but with filtered haplotype effects (LD filtered genotype)
outName='res_haplo_filtered'
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'
fi

# 4th order interaction scenarios
if [ "$order" = "4" ] ; then

# 2. baseline - epistatic truth (full genotype)
outName='res_epi'
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'

# 2b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_epi_PRS'
CheckIfNonLin_And_Lin $resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights'
fi

done  # end of numOrders  loop

done # end of numVersions  loop




####################
# Inference KNET:


order=1
b=0
phename='GWAS_'$order$'_v'$b
#  KNET, in the same 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 
phename='GWAS_'$order$'_v'$b

# additive scenarios
if [ "$order" = "1" ] ; then
echo "order is: "$order
# 1. baseline - additive truth (full genotype)
outName='res'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"


# 1b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_PRS'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --prs '$simsData$phename$knetPRSext$'  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"


# Haplotype effect tests (only order 1, ie additive)
# 3a. haplotype-effect - additive truth (incomplete coverage genotype)
outName='res_haplo'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST_haplo_'$b$' --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"


# 3b. same as a. but with PRS-CS premultiplied SNPs (incomplete coverage genotype)
outName='res_haplo_PRS'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST_haplo_'$b$' --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --prs '$simsData$phename$knetPRSext$'_haplo  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"

# 3c. same as a. but with filtered haplotype effects (LD filtered genotype)
outName='res_haplo_filtered'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST_haplo_filtered'$b$' --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"
fi

# 4th order interaction scenarios
if [ "$order" = "4" ] ; then
echo "order is: "$order

# 2. baseline - epistatic truth (full genotype)
outName='res_epi'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"


# 2b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_epi_PRS'
arguments='--out  '$resultsLoc$'sims/'$phename$'/'$outName$' knet --inference 1  --batch_size '$batch_size$' --gradient_batch_size '$gradient_batch_size$' --hyperopt 0 --earlystop 1 --oversampling 1 --cc 0 --gpu '$gpu$' --dropout '$dropout$' --firstLayerSize '$firstLayerSize$' --hidCount '$hidCount$' --hidAct '$hidAct$' --plink '$UKBB_PLINK1$'SIMS_TEST --pheno '$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno  --prs '$simsData$phename$knetPRSext$'  --loadWeights '$resultsLoc$'sims/'$phename$'/'$outName$'SaveWeights' 
#python3 /home/mk907/software/knet2/Knet.py $arguments
sbatch --cpus-per-task 3 --mem ${PLINKRAM}  --nodes=1 --gres=gpu:1 --time 12:0:0 --partition ampere --account DANESH-SL3-GPU -e $resultsLoc$'sims/'$phename$'/'$outName$'err' -o $resultsLoc$'sims/'$phename$'/'$outName$'out' --wrap "python3 /home/mk907/software/knet2/Knet.py $arguments"

fi

done  # end of numOrders  loop

done # end of numVersions  loop






# Check if all results are there...
#  KNET, in the same 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 
phename='GWAS_'$order$'_v'$b

# additive scenarios
if [ "$order" = "1" ] ; then
echo "order is: "$order
# 1. baseline - additive truth (full genotype)
outName='res'
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName

# 1b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_PRS'
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName

# Haplotype effect tests (only order 1, ie additive)
# 3a. haplotype-effect - additive truth (incomplete coverage genotype)
outName='res_haplo'
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName

# 3b. same as a. but with PRS-CS premultiplied SNPs (incomplete coverage genotype)
outName='res_haplo_PRS'
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName

# 3c. same as a. but with filtered haplotype effects (LD filtered genotype)
outName='res_haplo_filtered'
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName
fi

# 4th order interaction scenarios
if [ "$order" = "4" ] ; then
echo "order is: "$order
# 2. baseline - epistatic truth (full genotype)
outName='res_epi'
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName

# 2b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_epi_PRS'
Check_TEST_exists $resultsLoc$'sims/'$phename$'/'$outName

fi

done  # end of numOrders  loop

done # end of numVersions  loop

###############


b=0
order=1
phename='GWAS_'$order$'_v'$b
mkdir -p $resultsLoc$phename$'/eval/'
#  KNET, in the same 3 (+4) scenarios: 
for ((b=0; b<numSimReplicates; b++)); do 
for ((order=1; order<=numOrders; order++)); do 
phename='GWAS_'$order$'_v'$b
mkdir -p $resultsLoc$phename$'/eval/'



# additive scenarios
if [ "$order" = "1" ] ; then
echo "order is: "$order
# 1. baseline - additive truth (full genotype)
outName='res'
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'


# 1b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_PRS'
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'


# Haplotype effect tests (only order 1, ie additive)
# 3a. haplotype-effect - additive truth (incomplete coverage genotype)
outName='res_haplo'
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'


# 3b. same as a. but with PRS-CS premultiplied SNPs (incomplete coverage genotype)
outName='res_haplo_PRS'
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'

# 3c. same as a. but with filtered haplotype effects (LD filtered genotype)
outName='res_haplo_filtered'
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'
fi

# 4th order interaction scenarios
if [ "$order" = "4" ] ; then
echo "order is: "$order
# 2. baseline - epistatic truth (full genotype)
outName='res_epi'
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'

# 2b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_epi_PRS'
GenAndEval_ALL_PRS $resultsLoc$'sims/'$phename$'/'$outName $simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno' '0' $resultsLoc$phename$'/eval/'$outName '0'

fi

done  # end of numOrders  loop

done # end of numVersions  loop


#############

# extract all results into 3 text files ( we work off the correlation not correlation^2)
header="replicate\tres\tres_PRS\tres_haplo\tres_haplo_PRS\tres_haplo_filtered\tres_epi\tres_epi_PRS"
echo -e $header > $resultsLoc$'simRes_nonlinear'
echo -e $header > $resultsLoc$'simRes_noAct'
echo -e $header > $resultsLoc$'simRes_linear'

for ((b=0; b<numSimReplicates; b++)); do 

order=1
phename='GWAS_'$order$'_v'$b
# 1. baseline - additive truth (full genotype)
outName='res'
res_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# 1b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_PRS'
res_PRS_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_PRS_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_PRS_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# Haplotype effect tests (only order 1, ie additive)
# 3a. haplotype-effect - additive truth (incomplete coverage genotype)
outName='res_haplo'
res_haplo_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_haplo_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_haplo_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# 3b. same as a. but with PRS-CS premultiplied SNPs (incomplete coverage genotype)
outName='res_haplo_PRS'
res_haplo_PRS_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_haplo_PRS_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_haplo_PRS_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# 3c. same as a. but with filtered haplotype effects (LD filtered genotype)
outName='res_haplo_filtered'
res_haplo_filtered_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_haplo_filtered_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_haplo_filtered_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# 4th order interaction scenarios
order=4
phename='GWAS_'$order$'_v'$b

# 2. baseline - epistatic truth (full genotype)
outName='res_epi'
res_epi_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_epi_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_epi_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


# 2b. same as a. but with PRS-CS premultiplied SNPs (full genotype)
outName='res_epi_PRS'
res_epi_PRS_linear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_linear_r2')
res_epi_PRS_nonlinear=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_nonlinear_r2')
res_epi_PRS_noAct=$(awk '{ if(FNR == 2) {print $1 } }' $resultsLoc$phename$'/eval/'$outName$'_noAct_r2')


echo -e $b$"\t"$res_nonlinear$"\t"$res_PRS_nonlinear$"\t"$res_haplo_nonlinear$"\t"$res_haplo_PRS_nonlinear$"\t"$res_haplo_filtered_nonlinear$"\t"$res_epi_nonlinear$"\t"$res_epi_PRS_nonlinear >> $resultsLoc$'simRes_nonlinear'
echo -e $b$"\t"$res_linear$"\t"$res_PRS_linear$"\t"$res_haplo_linear$"\t"$res_haplo_PRS_linear$"\t"$res_haplo_filtered_linear$"\t"$res_epi_linear$"\t"$res_epi_PRS_linear >> $resultsLoc$'simRes_linear'
echo -e $b$"\t"$res_noAct$"\t"$res_PRS_noAct$"\t"$res_haplo_noAct$"\t"$res_haplo_PRS_noAct$"\t"$res_haplo_filtered_noAct$"\t"$res_epi_noAct$"\t"$res_epi_PRS_noAct >> $resultsLoc$'simRes_noAct'

done # end of numVersions  loop


# !!! HERE
#############################################
# VI. Plot Results
#############################################


#$resultsLoc$'simRes_nonlinear'
#$resultsLoc$'simRes_noAct'
#$resultsLoc$'simRes_linear'

arguments='/home/mk907/scripts/R/PRS_GXE_SimsPlotter.R '$resultsLoc$'simRes '$resultsLoc$'simRes_res PRS_GxE_simulations'
Rscript $arguments
head $resultsLoc$'simRes_res_tab'
wc -l $resultsLoc$'sims/GWAS_1_res_haplo_r2' # this has the 20 PRS r2s
wc -l $resultsLoc$'sims/GWAS_1_res_haplo_r2_redux'  # this has the 20 PRS r AND r2s too

##results: 
# scenario	            mean_nonlinear	    t-test_p
# res_epi	            0.108951012868357	7.27079503864627e-05
# res	                0.0366589894767511	4.37273552878985e-10
# res_haplo	            0.0163957566573847	0.000272450943566516
# res_haplo_PRS	        0.00238875916613838	0.398304785133278
# res_haplo_filtered   -0.0121123686785715	4.40292384280496e-11
# res_epi_PRS	        0.0947281379734188	1.74860099845811e-06

# a significant positive value means that the non-linear model is better

# res_epi: positive, the expected, when there is true epistasis the NN will find it (but this is still subject to haplotype effects, as we had 0.03665, for the additive scenario WITHOUT haplotype effects) 
# res_epi_PRS: positive, and somewhat smaller than res_epi, this is expected, as it seems that at least part of res_epi was due to haplotype effects! 

# res: positive, NN will find non-linearity, even when no true non-linearity exists, this is somewhat unexpected, as we have complete coverage here, IE no haplotype effects
# res_haplo: positive,  haplotype effects: NN finds non-linearity, due to haplotype effects (but this is lower than 'res', probably due to lower power, as we have removed half the causal variants, so not entirely unexpected.)
# res_haplo_PRS: zero (non-significant), seems that per-SNP PRS correction eliminates haplotype-effects
# res_haplo_filtered: negative, ie linear model is better, which suggests that this strategy overcorrects against haplotype effects (as this should also be 0)

# scenario	        mean_nonlinear	    t-test_p
# res_PRS_adv	    0.651112236240643	9.11211355917109e-19
# res_epi_PRS_adv	0.958950762878411	7.01357510815655e-19
# So the '_PRS' effect always improves upon the baseline in absolute terms, which suggests that this is an advantageous strategy



# difference between a,b and c tell if/how we can mitigate the haplotype effect artifact




# 1a. baseline - additive truth (full genotype)                                     res_nonlinear
# 1b. same as a. but with PRS-CS premultiplied SNPs (full genotype)                 res_PRS_nonlinear
# 2a. baseline - epistatic truth (full genotype)                                    res_epi_nonlinear
# 2b. same as a. but with PRS-CS premultiplied SNPs (full genotype)                 res_epi_PRS_nonlinear
# 3a. haplotype-effect - additive truth (incomplete coverage genotype)              res_haplo_nonlinear
# 3b. same as a. but with PRS-CS premultiplied SNPs (incomplete coverage genotype)  res_haplo_PRS_nonlinear
# 3c. same as a. but with filtered haplotype effects (LD filtered genotype)         res_haplo_filtered_nonlinear




# Plots: 

# A) NON-LINEARITY: (using 2a, res_epi_nonlinear): if the nonlinearity works this shows as a proof of concept that NNs CAN work, if there is non-linearity
# B) HAPLO-TYPE EFFECTS: (using 1a and 3a/b/c):
# 1a, res_nonlinear: here the nonlinear model should NOT do better than the linear model, as we have full coverage
# 3a, res_haplo_nonlinear: here the nonlinear method SHOULD do better, even though we have an additive truth, due to the haplotype effects and incomplete coverage
# this shows the problems that naive NN PRS projects have had
# 3b, res_haplo_PRS_nonlinear. the non-linear method should hopefully not do better
# this shows if haplotype effects could be mitigated by using per-SNP PRS weights (not because PRS weights gain any non-linear capability, but because they model LD and thereby already 'impute'
# (this should also be better than 3a)
# 3c, res_haplo_filtered_nonlinear. the non-linear method should hopefully not do better
# this shows if haplotype effects could be mitigated by using my LD/distance filtering method
# this will be WORSE than 3a, as we will only have very few SNPs, but at least get a 'clean' result

# so there will be 5 independent boxplots / table results 


# genotypes
$UKBB_PLINK1$'SIMS_TRAINVALID'
$UKBB_PLINK1$'SIMS_TEST'

$UKBB_PLINK1$'SIMS_TRAINVALID_haplo_'$b
$UKBB_PLINK1$'SIMS_TEST_haplo_'$b


$UKBB_PLINK1$'SIMS_TRAINVALID_haplo_filtered'$b
$UKBB_PLINK1$'SIMS_TEST_haplo_filtered'$b


# phenotypes
$simsDataTruthset$'phenos/GWAS_KNET_trainvalid_'$order$'_v'$b$'.pheno'
$simsDataTruthset$'phenos/GWAS_KNET_test_'$order$'_v'$b$'.pheno'


# additive PRS file for the full data
$simsData$'PRSProfiles/'$phename$'_PRSCS_all.sscore'


# additive PRS file for the haplotype effect data
$simsData$'PRSProfiles/'$phename$'_haplo_PRSCS_all.sscore'

$simsData$'PRSProfiles/'$phename$'_haplo_PRSCS_all.sscore'