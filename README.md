
# Code repository for "Performance of deep-learning based approaches to improve polygenic scores"

DOI: 10.1101/2024.10.23.24315973

This respository represents the last snapshot of the bash and R scripts used to generate our results and is provided as-is. As the analysis involved a lot of input/output operations on very large files, these were generated asynchronously on a cluster. Thus these scripts are meant to be executed on the command line manually, block-by-block, waiting for the remote jobs to finish and verifying the integrity of the resulting files at each step. To reduce code duplication, certain functions that were reused multiple times are defined only once across all files, however, they may be called from different scripts.

**The simulation analyses:**
1. PRS_GXE_sims.sh: code relevant for the simulation analyses

**The real data analyses:**
1. PRS_GXE_v2.sh: all real data analyses
2. PRS_GX_functions_v2.sh: common functions

**neural-network pytorch model:**
1. scripts/python/Knet/: Visual studio project that implements the neural-network models used in the study

**revisions:**
1. scripts under subfolders /rev1: contain the scripts for the new analyses for the first revision.

**REQUIREMENTS & INSTALLATION:**
1. Install Pytorch (v1.9.0+cu111): https://pytorch.org/get-started/
(for installation times see the Pytorch documentation)
2. Copy the contents of /scripts/python/Knet/ to somewhere conventient

**DEMO:**

The repository inludes some toy data as a minimal example that should take a few seconds to run:

*yourDir='<YOURPATH>/demo/'*

*knet='<YOURPATH>/Knet.py'*

**TRAINING**

*python3 $knet --out $yourDir$'results/' knet --batch_size 4 --gradient_batch_size 4 --hyperopt 0 --cc 0 --gpu 0 --firstLayerSize 10 --hidCount 1 --hidAct 4 --plink $yourDir$'/trainvalid' --pheno $yourDir$'/trainValid.pheno' --validSet $yourDir$'/valid' --saveWeights $yourDir$'results/SaveWeights'*

generates and saves model weights with prefix 'SaveWeights' and diagnostic results for the trainValid set 

**Generate PRS for the test set**

*python3 $knet --out $yourDir$'results/NN_' knet --inference 1 --batch_size 4 --gradient_batch_size 4 --hyperopt 0 --cc 0 --gpu 1 --firstLayerSize 10 --hidCount 1 --hidAct 4 --plink $yourDir$'/test' --pheno $yourDir$'/test.pheno' --loadWeights $yourDir$'results/SaveWeights'*

**3 PRS outputs:**

**NN_yhat_TEST_noAct_retrain.txt:** the linear NN PRS for the test set indis

**NN_yhat_TEST.txt:** the non-linear NN PRS for the test set indis

**NN_yhat_TEST_noAct.txt:** alternative version the linear NN PRS for the test set indis (this is not practical as it just switches off the activation without retraining the model)

and

**NN_FIDs_TEST.txt:** The IDs of your indis in the same order as the PRS
