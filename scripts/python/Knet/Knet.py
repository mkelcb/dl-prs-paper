# -*- coding: utf-8 -*-

#MIT License

#Copyright (c) 2017 Marton Kelemen

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.





#from com.application.logic.knet import knet_manager
import knet_manager_pytorch
#from com.application.logic.knet import knet_manager_keras

# https://docs.python.org/3/library/argparse.html
#https://docs.python.org/3/howto/argparse.html 
import argparse
import os
import gc


def set_Threads(args) : 
    if args.threads is not None  :
        os.environ['MKL_NUM_THREADS'] = args.threads # '16'  # use here the N, where N: the number of cores acailable or limit to 1
        os.environ['MKL_DYNAMIC'] = 'FALSE'
        os.environ['OMP_NUM_THREADS'] = '1'
        
        print("set MKL number of threads to: " + str(args.threads))
   
    
def set_nixMem(args) : 
    if args.nixMem is not None  :
        import resource # this only exists on Unix/Linux based systems
        rsrc = resource.RLIMIT_AS
        soft, hard = resource.getrlimit(rsrc)
        print('Soft limit starts as  :', soft)
        print('Hard limit starts as  :', hard)
        
        resource.setrlimit(rsrc, (args.nixMem * 1048576, hard)) #limit
        
        soft, hard = resource.getrlimit(rsrc)
        print('Soft limit changed to :', soft)
        print('Hard limit changed to  :', hard)


def runKnet(args) :
    print('Knet Neural net started')
    set_Threads(args)
    set_nixMem(args) 
    args.earlystop = args.earlystop == 1
    args.bnorm = args.bnorm == 1
    args.half = args.half == 1
    knet_manager_pytorch.runKnet(args)  
    gc.collect() 


##################################################################################################
# setup COmmand line parser
##################################################################################################

parser = argparse.ArgumentParser()

# overall
parser.add_argument("--out",required=True, help='an output location is always required')
parser.add_argument("--threads",required=False, help='set number of threads used by multithreaded operations')
parser.add_argument("--nixMem",required=False, type=int, help='Memory limit for *nix based systems in Megabytes')

subparsers = parser.add_subparsers()
subparsers.required = True
subparsers.dest = 'either knet, scanner, h2, kinship, kinmerge or merge' # hack to make subparser required


# create the parser for the "a" command
parser_knet = subparsers.add_parser('knet')
parser_knet.add_argument('--plink', help='A plink genotype file.') # the location of the train set binaries
parser_knet.add_argument("--pheno", required=True)
parser_knet.add_argument("--prs", help='A polygenic score file that will be used to weight the SNPs. It should have structure (with header): hm_chr	hm_pos	effect_allele	other_allele	effect_weight')
parser_knet.add_argument("--device",required=False, help='the GPU device used to host the master copy of the model, default 0', default=0, type=int)   


parser_knet.set_defaults(func=runKnet)

# knet subparams
parser_knet.add_argument("--covars_IDs", help='The Individual IDs for the covariates')
parser_knet.add_argument("--covars_cont", help='Continuous covariates. Headerless file, where the columns are aligned to covars_IDs') 
parser_knet.add_argument("--covars_factor", help='Factor covariates. Headerless file, where the columns are aligned to covars_IDs') 

parser_knet.add_argument("--prs_indi", help='(optional) individual level PRS file with signature: IID\tPHENO1\tSCORE1_SUM') 

parser_knet.add_argument("--validSet") # the location for list of samples which are meant to be set aside for validation

parser_knet.add_argument("--loadWeights") # from where we want to load the weights
parser_knet.add_argument("--saveWeights") # where we wnt to save weights
parser_knet.add_argument("--savFreq", default=-1, type=int) # how frequently we make backups of Weights
parser_knet.add_argument("--epochs", default=20, type=int) # how many epochs
parser_knet.add_argument("--momentum",required=False, help='momentum used for the optimizer. default is 0.9', default=0.9, type=float)        # 
parser_knet.add_argument("--learnRate",required=False, help='learnRate used for the optimizer. default is 0.001', default=0.001, type=float)        # 
parser_knet.add_argument("--LRdecay",required=False, help='Learning rate decay, default 0.96 (to disable set it to -1)', default=-1, type=float)        # 


#parser_knet.add_argument("--validPhen") # the location for the binaries for the validation set phenotypes
                  
parser_knet.add_argument("--cc", type=int)  # ,required=False  # if phenotype is case control
parser_knet.add_argument("--recodecc", type=int)  # ,required=False       # if we want to recode case control to quantitative
parser_knet.add_argument("--randomSeed", default=1, type=int)                       
parser_knet.add_argument("--hidCount", default=0, type=int)     # number of hidden layers
  
parser_knet.add_argument("--hidAct", default=0, type=int, help=' the hidden layer activations ( 1 = sigmoid, 2 = RELU, 3 = linear, 4 = softplus, 5 =  LeakyReLU, 6 =SELU)')        #
         


parser_knet.add_argument("--gradient_batch_size",required=False, help='effective size of minibatches used for gradient calculation, default :32', default=32, type=int)        # 
    
parser_knet.add_argument("--batch_size", default=32, type=int, help='the size of the minibatches, default :32')          # the size of the minibatches, use 0 for no minibatches (IE train all at once)
parser_knet.add_argument("--bnorm",required=False, help='if batchnorm (1, default) or group norm is to be used ', default=1, type=int)   
parser_knet.add_argument("--lr_decay", default=0, type=int)        # learning rate decay should be enabled (1) or not (0)
parser_knet.add_argument("--optimizer", default=0, help='the optimizer, 0 for SGD (the default), 1 for ADAM, and 2 for AMSGrad', type=int)        # 
#parser_knet.add_argument("--float", default=32, type=int)        # the float precision, valid options are 16, 32 and 64 (the default)
parser_knet.add_argument("--inference", default=0, type=int)        # if an Inference (ie deep dreaming) run is to be performed (1) or training run should be performed (0)
parser_knet.add_argument("--orig", default=0, type=int)
parser_knet.add_argument("--firstLayerSize", default=100, type=int) # the number of units in the first layer
parser_knet.add_argument("--dropout", default=-1, type=float)  # % of units to switch off at each iteration

#parser_knet.add_argument("--snpIndices") 
#parser_knet.add_argument("--mns") 
#parser_knet.add_argument("--sstd") 
#parser_knet.add_argument("--snpIDs") 

parser_knet.add_argument("--convLayers", default=0, type=int, help='how many convolutional layers to add (0 for disabled)') # 
parser_knet.add_argument("--convFilters", default=500, type=int, help ='the number of filters that we will use in the first layer, each subsequent layer will have i * this many filters ') # 
parser_knet.add_argument("--widthReductionRate", default=1, help='The rate at which the network "thins" IE if we start at 1000 neurons in layer 1, then at rate of 1 (default), we half it every layer, with a rate of 2, it will half every second layer Ie we will get two layers with 1000 units each, and then two 500 units etc', type=int) 
parser_knet.add_argument("--half",required=False, help='if FP16 should be used (default no)', default=0, type=int)   
parser_knet.add_argument("--gpu",required=False, help='the number of gpus to be used. 0 for cpu.', default=0, type=int)              
                        
# parser_knet.add_argument("--topology", required=True) # the location of the file that describes the network's topology (IE number and size of layers etc)
parser_knet.add_argument("--predictPheno", default=-1, type=int) # if network should save phenotype predictions to a location at the end, for a validation set                              
parser_knet.add_argument("--num_CPU", default=1, type=int) # the number of CPU cores Keras should use                  
parser_knet.add_argument("--qc", default=0, type=int, help='if SNP QC is to be performed (1) or no (0)') # 

parser_knet.add_argument("--l2", default=0.0, type=float)        # the L2 regularizer shrinkage param     
    
parser_knet.add_argument("--hyperopt", default=0, type=int, help='if best parameter settings are to be found via hyperopt semi-random search, 0 for NO hyperopt, otherwise the number of trials') # 

parser_knet.add_argument("--epochMaxImproveThreshold",required=False, help='Max number of epochs until no improvement before stopping. default is 12', default=12, type=int)        #  
parser_knet.add_argument("--earlystop",required=False, help='if early stop mechanism is to be applied (default True)', default=1, type=int)

parser_knet.add_argument("--linearInference", default=0, type=int) # if the activation functions should be switched off (1) for inference or not (0)
parser_knet.add_argument("--oversampling", help ='If oversampling logic is enabled for when there is an imbalance between cases and controls') 

parser_knet.add_argument("--disablePRS", help ='If the PRS weights are disabled. (useful to make sure we use same SNPs as if had a PRS)') 

parser_knet.add_argument("--redoEarlyConv", default=1, type=int, help='if models that converged at epochs 0 or 1 are redone (1, the default) or not (0)') # 


# retreive command line arguments
args = parser.parse_args()
args.func(args)

# toy test
# python /nfs/users/nfs_m/mk23/software/knet/knet.py --out /nfs/users/nfs_m/mk23/test/pytest/toyregions_ --threads 2 scanner --scanner /nfs/users/nfs_m/mk23/data/gwas2/toy/wtccc2_hg19_toy --pheno /nfs/users/nfs_m/mk23/data/gwas2/toy/wtccc2_hg19_toy.pheno --saveEigSum /nfs/users/nfs_m/mk23/test/pytest/toyeig



##################################################################################################
##################################################################################################
# Local Tests



## KNET MAIN
# args = parser.parse_args(['--out', 'C:/0LocalHPC/PRS_GXE/dummy/out/dummy_out','knet', '--plink','C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK', '--pheno', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_bin.pheno', '--cc', '1', '--covars_cont', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_numericBinary', '--covars_factor', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_factors', '--covars_IDs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_indiIDs', '--prs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS',  '--validSet', 'C:/0LocalHPC/PRS_GXE/dummy/validIndis',  '--batch_size', '32' ,  '--hyperopt', '7' ,  '--earlystop', '1', '--saveWeights', 'C:/0LocalHPC/PRS_GXE/dummy/out/SaveWeights' ]) # -- --  = 32  , '--loadEigSum', '../../../0cluster/data/data/toy/eig/'

# 
# args = parser.parse_args(['--out', 'C:/0LocalHPC/PRS_GXE/dummy/out/dummy_out','knet', '--plink','C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK', '--pheno', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_missing.pheno', '--covars_cont', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_numericBinary', '--covars_factor', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_factors', '--covars_IDs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_indiIDs', '--prs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS',  '--validSet', 'C:/0LocalHPC/PRS_GXE/dummy/validIndis',  '--batch_size', '32' ,  '--hyperopt', '7' ,  '--earlystop', '1', '--saveWeights', 'C:/0LocalHPC/PRS_GXE/dummy/out/SaveWeights' ]) # -- --  = 32  , '--loadEigSum', '../../../0cluster/data/data/toy/eig/'
# args = parser.parse_args(['--out', 'C:/0LocalHPC/PRS_GXE/dummy/out/dummy_out','knet', '--plink','C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK', '--pheno', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_missing.pheno', '--covars_cont', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_numericBinary', '--prs_indi', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PRS.sscore', '--covars_factor', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_factors', '--covars_IDs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_indiIDs', '--prs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS',  '--validSet', 'C:/0LocalHPC/PRS_GXE/dummy/validIndis',  '--batch_size', '32' ,  '--hyperopt', '7' ,  '--earlystop', '1', '--saveWeights', 'C:/0LocalHPC/PRS_GXE/dummy/out/SaveWeights' ]) # -- --  = 32  , '--loadEigSum', '../../../0cluster/data/data/toy/eig/'


# no plink
#args = parser.parse_args(['--out', 'C:/0LocalHPC/PRS_GXE/dummy/out/dummy_out','knet', '--pheno', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_missing.pheno', '--covars_cont', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_numericBinary', '--covars_factor', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_COVS_factors', '--covars_IDs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_indiIDs', '--prs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS',  '--validSet', 'C:/0LocalHPC/PRS_GXE/dummy/validIndis',  '--batch_size', '32' ,  '--hyperopt', '7' ,  '--earlystop', '1', '--saveWeights', 'C:/0LocalHPC/PRS_GXE/dummy/out/SaveWeights' ]) # -- --  = 32  , '--loadEigSum', '../../../0cluster/data/data/toy/eig/'
# no covars
#args = parser.parse_args(['--out', 'C:/0LocalHPC/PRS_GXE/dummy/out/dummy_out','knet', '--plink','C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK', '--pheno', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PLINK_missing.pheno', '--prs', 'C:/0LocalHPC/PRS_GXE/dummy/dummyPRS',  '--validSet', 'C:/0LocalHPC/PRS_GXE/dummy/validIndis',  '--batch_size', '32' ,  '--hyperopt', '7' ,  '--earlystop', '1', '--saveWeights', 'C:/0LocalHPC/PRS_GXE/dummy/out/SaveWeights' ]) # -- --  = 32  , '--loadEigSum', '../../../0cluster/data/data/toy/eig/'



