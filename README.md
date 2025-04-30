
# Code repository for "Performance of deep-learning based approaches to improve polygenic scores"

DOI: 10.1101/2024.10.23.24315973

This respository represents the last snapshot of the bash, R and python scripts used to generate our results and is provided as-is. As the analysis involved a lot of input/output operations on very large files, these were generated asynchronously on a cluster. Thus these scripts are meant to be executed on the command line manually, block-by-block, waiting for the remote jobs to finish and verifying the integrity of the resulting files at each step. To reduce code duplication, certain functions that were reused multiple times are defined only once across all files, however, they may be called from different scripts.

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
2. Install rest of dependencies listed in the 'requirements.txt', which can be done by Visual Studio:
https://learn.microsoft.com/en-us/visualstudio/python/managing-required-packages-with-requirements-txt
3. Copy the contents of /scripts/python/Knet/ to somewhere conventient

**DEMO:**

The repository inludes some toy data as a minimal example that should take a few seconds to run:

*yourDir='<YOURPATH>/demo/'*

*knet='<YOURPATH>/Knet.py'*

*# (knet='scripts/python/Knet/Knet.py' when running from within the 'dl-prs-paper-master' folder.)*

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

**Detailed help for all parameters can be displayed by:**

*python3 $knet knet --help*

<code>  -h, --help            show this help message and exit
  --plink PLINK         A plink genotype file.
  --pheno PHENO
  --prs PRS             A polygenic score file that will be used to weight the
                        SNPs. It should have structure (with header): hm_chr
                        hm_pos effect_allele other_allele effect_weight
  --device DEVICE       the GPU device used to host the master copy of the
                        model, default 0
  --covars_IDs COVARS_IDS
                        The Individual IDs for the covariates
  --covars_cont COVARS_CONT
                        Continuous covariates. Headerless file, where the
                        columns are aligned to covars_IDs
  --covars_factor COVARS_FACTOR
                        Factor covariates. Headerless file, where the columns
                        are aligned to covars_IDs
  --prs_indi PRS_INDI   (optional) individual level PRS file with signature:
                        IID PHENO1 SCORE1_SUM
  --validSet VALIDSET
  --loadWeights LOADWEIGHTS
  --saveWeights SAVEWEIGHTS
  --savFreq SAVFREQ
  --epochs EPOCHS
  --momentum MOMENTUM   momentum used for the optimizer. default is 0.9
  --learnRate LEARNRATE
                        learnRate used for the optimizer. default is 0.001
  --LRdecay LRDECAY     Learning rate decay, default 0.96 (to disable set it
                        to -1)
  --cc CC
  --recodecc RECODECC
  --randomSeed RANDOMSEED
  --hidCount HIDCOUNT
  --hidAct HIDACT       the hidden layer activations ( 1 = sigmoid, 2 = RELU,
                        3 = linear, 4 = softplus, 5 = LeakyReLU, 6 =SELU)
  --gradient_batch_size GRADIENT_BATCH_SIZE
                        effective size of minibatches used for gradient
                        calculation, default :32
  --batch_size BATCH_SIZE
                        the size of the minibatches, default :32
  --bnorm BNORM         if batchnorm (1, default) or group norm is to be used
  --lr_decay LR_DECAY
  --optimizer OPTIMIZER
                        the optimizer, 0 for SGD (the default), 1 for ADAM,
                        and 2 for AMSGrad
  --inference INFERENCE
  --orig ORIG
  --firstLayerSize FIRSTLAYERSIZE
  --dropout DROPOUT
  --convLayers CONVLAYERS
                        how many convolutional layers to add (0 for disabled)
  --convFilters CONVFILTERS
                        the number of filters that we will use in the first
                        layer, each subsequent layer will have i * this many
                        filters
  --widthReductionRate WIDTHREDUCTIONRATE
                        The rate at which the network "thins" IE if we start
                        at 1000 neurons in layer 1, then at rate of 1
                        (default), we half it every layer, with a rate of 2,
                        it will half every second layer Ie we will get two
                        layers with 1000 units each, and then two 500 units
                        etc
  --half HALF           if FP16 should be used (default no)
  --gpu GPU             the number of gpus to be used. 0 for cpu.
  --predictPheno PREDICTPHENO
  --num_CPU NUM_CPU
  --qc QC               if SNP QC is to be performed (1) or no (0)
  --l2 L2
  --hyperopt HYPEROPT   if best parameter settings are to be found via
                        hyperopt semi-random search, 0 for NO hyperopt,
                        otherwise the number of trials
  --epochMaxImproveThreshold EPOCHMAXIMPROVETHRESHOLD
                        Max number of epochs until no improvement before
                        stopping. default is 12
  --earlystop EARLYSTOP
                        if early stop mechanism is to be applied (default
                        True)
  --linearInference LINEARINFERENCE
  --oversampling OVERSAMPLING
                        If oversampling logic is enabled for when there is an
                        imbalance between cases and controls
  --disablePRS DISABLEPRS
                        If the PRS weights are disabled. (useful to make sure
                        we use same SNPs as if had a PRS)
  --redoEarlyConv REDOEARLYCONV
                        if models that converged at epochs 0 or 1 are redone
                        (1, the default) or not (0)
</code>