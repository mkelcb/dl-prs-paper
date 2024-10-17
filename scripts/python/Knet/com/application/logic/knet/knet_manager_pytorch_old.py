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

import gc
from math import remainder
import numpy as np
import pandas as pd
from numpy.linalg import norm 
from scipy import stats
from pathlib import Path
import random
import os
import time
import sys 
import matplotlib
import platform
if platform.system().find('Windows') == -1 :
    print("Matlab uses Agg as we are on a *nix")
    matplotlib.use('Agg')

import matplotlib.pyplot as plt   
from functools import partial
from types import SimpleNamespace
import copy
from ast import literal_eval
from random import choices
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.nn.init as init
from torch.autograd import Variable

from sklearn.model_selection import train_test_split
#from sklearn.metrics import accuracy_score, cohen_kappa_score, roc_auc_scores
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials, space_eval
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

        
from ....application.utils.plotgen import exportNNPlot
from ....application.utils.geno_qc import removeList, genoQC_all, standardise_Genotypes, getSizeInMBs , zscore
from ....application.logic.knet.knet_main_pytorch import weight_init, EPSILON, learn, setModelMode , NETWORK_DATATYPE, getNetworkDatatype_numpy, getModel, isLayerActivation, getMiniBatchData,loadEpochStats
from ....io.knet_IO import loadPLINKPheno, loadPLINK, loadMatrixFromDisk, loadVectorFromDisk, get_d_code, loadNumericFromTextFile, loadTextList, matchIndices, loadPRSFile, matchIndices2, validTrainSplit, writeVectorToDisk
from ....io import knet_IO
#from ....application.utils.sumstatsProcessor import loadRegionsData


#from ....application.logic.knet.knet_association import mainAssoc

# delta = (Ve/Vg)
# delta = (1-h2) / h2
#args, args.epochs, args.learnRate, args.momentum, args.evalFreq, args.savFreq, args.predictPheno, args.loadWeights, args.saveWeights, args.randomSeed, args.hidCount, args.hidl2, args.hidAct

#args = parser.parse_args(['--out', 'C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/tests/0pytorch_tests/' ,'knet', '--knet', 'C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/data/genetic/short', '--pheno', 'C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/data/genetic/phenew.pheno.phen',  '--epochs', '21', '--learnRate', '0.00005', '--momentum', '0.9', '--evalFreq', '1',  '--recodecc'    , '0' ,  '--hidCount'    , '5' ,  '--hidl2'    , '1.0' ,  '--hidAct'    , '2' , '--cc', '0' ,'--inference', '0'   ]) # ,'--loadWeights', 'C:/0Datasets/NNs/genetic/weights/' ,'--snpIndices', 'C:/0Datasets/NNs/genetic/nn_SNPs_indices.txt' ,'--mns', 'C:/0Datasets/NNs/genetic/data_mns','--sstd', 'C:/0Datasets/NNs/genetic/data_sstd','--snpIDs', 'C:/0Datasets/NNs/genetic/nn_SNPs.txt'

    
###############################################################################
# Global vars
###############################################################################
device = None

y= None
y_valid= None
pheno_is_binary = False

best_epoch_string =  "best_epoch.txt"
_results =  "_results"
_optimizer =  "_optimizer"
_stats =  "_stats"

###############################################################################
# Model construction
#region Model construction
###############################################################################
def build_model(args, numSNPs,  numCovars, num_y_classes, suppressPrint = False ) :
    torch.cuda.empty_cache()
    suppressPrint = False
    hLayerCount = args.hidCount
    BNEnabled = int(args.bnorm) == 1
    torch.manual_seed(args.randomSeed)
    layers = []
    lastLayerSize = args.firstLayerSize
    lastOutput = numSNPs + numCovars # + numCovarsFactor, the input size is the number of SNPs plus all the continuous covars
    

    for i in range(1,hLayerCount+1) : # iterate 1 based, otherwise we will get a reduction after the first layer, no matter the widthReductionRate, as 0 is divisible by anything
        layers.append( nn.Linear(lastOutput, lastLayerSize) ) # torch layers are parametrized as: In , Out (IE num rows, num cols (IE the number of neurons))
        
        if BNEnabled : layers.append( nn.BatchNorm1d(lastLayerSize ) )   
        addActivation(layers,args.hidAct)
        if args.dropout != -1 : addDropout(layers,args)
        
        if suppressPrint == False : print("added layer at depth: " + str(i) + " with width: " + str(lastLayerSize))
        
        lastOutput = lastLayerSize
        # control the 'fatness' of the network: we reduce the width at a given rate: if this is 1, then at every subsequent layer, if its 2, then every 2nd layer etc
        if i % args.widthReductionRate == 0 :  lastLayerSize = lastLayerSize // 2
        if lastLayerSize < 2 : break # if 
        
        
    if num_y_classes > 1 : 
        layers.append( nn.Linear(lastOutput, num_y_classes) )
        layers.append( nn.Softmax() )
    else :
        layers.append( nn.Linear(lastOutput, 1) )

    # initialize the model: https://stackoverflow.com/questions/49433936/how-to-initialize-weights-in-pytorch
    model = nn.Sequential(*layers)
    model.apply(weight_init)

    
    return(model)
      #np.sqrt(0.5) * 0.3

def addActivation(layers, hidAct): 
    if hidAct == 1 : layers.append( nn.Sigmoid() )	
    elif hidAct == 2 : layers.append( nn.ReLU()	)
    elif hidAct == 5 : layers.append( nn.LeakyReLU(negative_slope=0.001) )	
    elif hidAct == 4 : layers.append( nn.Softplus()	)
    elif hidAct == 6 : layers.append( nn.SELU()	)
   # elif hidAct == 3 :  print("no activatioN")


class Flatten(nn.Module):
    def forward(self, input):
        #print("input hape is:", input.shape, flush=True)
        return input.view(input.size(0), -1)

    
def addDropout(layers, args): 
    if args.hidAct == 6 : layers.append( nn.AlphaDropout(p=args.dropout)	) # use SELU for alpha dropout
    else : layers.append( nn.Dropout(p=args.dropout)	)
    
 
def model_diagnostics(model, args):
    totalParams = sum([p.numel() for p in model.parameters()])
    totalParams # 726270  # 211
    resultsLines = []
    resultsLines.append("model has total num params: " + str( totalParams) )

    # a better view is to aggregate the number of model params per the group of hidden layers, IE the linear, BN, Dropout and activation together
    fcLayerParams = [] 
    mayorLayerCounter = 0
    for i in range(len( model)) :
        paramsPerLayer = 0
        currentLayer = model[i]
        if type(currentLayer) == nn.Linear and i != 0 or i == len( model) -1:  # we dont keep track of 'major' layers in this, so we just infer them by checking if we encountered a new 'linear' layer
            if i == len( model) -1 : paramsPerLayer += sum([p.numel() for p in currentLayer.parameters()]) # if it is the last layer we want to add get its parameters
            
            fcLayerParams.append(paramsPerLayer)
            resultsLines.append("FC layer "  +str(mayorLayerCounter) + " has num params: " +str(paramsPerLayer) + " / "  +str( np.round( (paramsPerLayer/ totalParams) * 100, 2) ) + "%" )
            fcLayerParams.append(paramsPerLayer)
            mayorLayerCounter += 1
            paramsPerLayer = 0
            
            
        if i != len( model) -1 : paramsPerLayer += sum([p.numel() for p in currentLayer.parameters()])

    
    allFCParamsFound = np.sum(fcLayerParams  ) 
    totalFound = allFCParamsFound
 
    resultsLines.append("accounted for all params: " + str( totalParams == totalFound) )
    with open(args.out+"_diag", "w") as file: 
        for i in range(len(resultsLines) ):
            print(resultsLines[i])
            file.write(str(resultsLines[i])  + "\n")
            
    print("written diagnostics info to", args.out+"_diag")  
#endregion
    
###############################################################################
# Utils
#region Utils
###############################################################################   
    
def dataLoad(dataLoc, loadPhenos = True, replaceMissing = False) : # generic input data loading function that works for both PLINK binaries as well as LDPred PRS binaries
    print("KNET: loading plink file")
    genotypeData = loadPLINK(dataLoc, loadPhenos = loadPhenos, replaceMissing = replaceMissing) 
    genotypeData["plink"] = True
    return(genotypeData)


def modelAlreadyRun (args, postFix = "") : # util, to determine if a model has already run,
    if args.saveWeights is not None :    
        pastModel = Path(args.saveWeights + postFix)
        if pastModel.is_file():
            print(postFix,"model already run")
            return(True)
        else : return(False)
    else : return(False)
 
    #suppliedWeights= training_weights
    #indices = remainderBatch
def meanLevelLoss(y, indices, num_batches =1) : #calculates the loss for a mean level prediction
    
    y_used = y[indices]
    mean_y = np.mean(y_used, axis=0, keepdims=True)
    meanPrediction = np.ones(y_used.shape, dtype=np.float32) * mean_y

  
    meanPrediction = torch.Tensor( meanPrediction)
    if pheno_is_binary : criterion = nn.CrossEntropyLoss() 
    else :  criterion = nn.MSELoss()
    
    y_torch =  torch.Tensor(y_used)
    
    print("meanPrediction is: ", meanPrediction.shape, flush=True) 
    print("y_torch is: ", y_torch.shape, flush=True) 


    
    loss =  criterion(meanPrediction,y_torch )
    # the loss is averaged along both axes for torch, so if we want the un-averaged version, we multiply by their product
    loss = loss.item()  * ( y_used.shape[0] *  y_used.shape[1] )  /num_batches  # calculate error
    return(loss)
    #loss_torch = criterion(meanPrediction,y_torch ).item() # 0.0012115546269342303
    #manualLoss = np.sum( (y -meanPrediction)**2 ) / ( y.shape[0] *  y.shape[1] ) # 0.001211554640314582



#endregion

###############################################################################
#region QC and Training
############################################################################### 
multiGPU = False     
def runKnet(args) :
    ###########################################################################
    # (I) LOAD: load, preprocess and minibatch the data
    # region (I) LOAD: load, preprocess and minibatch the data
    ##########
    linearPostFix = "_linear"
    suppressPrint = False
    global multiGPU; global pheno_is_binary 
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    y_mean = None ;  y_sd = None
    print("KNeT via Pytorch backend, version: ",torch.__version__)
    # default QC settings
    _minObserved = 0.95
    _minMAF = 0.01
    _minVariance = 0.02
    _standardiseData = True
    _replaceMissing = True

    # load plink binary / phenotypes want to load them here, so that we can release the memory once the raw data is no longer used
    cc = True
    if args.cc == 0 or args.cc is None : cc = False
    pheno_is_binary = cc

    recodecc = True
    if args.recodecc == 0 or args.recodecc is None: recodecc = False
    start = time.time()

    y, indisWithPhenos = loadPLINKPheno(args.pheno, caseControl = cc, recodeCaseControl = recodecc) 
    y = y.reshape(-1,1) # enforce 2D , otherwise z-score wont work

    ##################
    genotypeData = dataLoad(args.knet, loadPhenos = False) # knet_IO.loadPLINK(args.knet, loadPhenos = False) 
    #plinkVersion = genotypeData["plink"]

    M_orig = genotypeData["M"]
    irsIds = genotypeData["rsid"]
    IDs = genotypeData["IDs"] 
    IDs[0] = np.array(IDs[0]) ; IDs[1] = np.array(IDs[1])
    indicesKept = np.asarray( range(M_orig.shape[1]) )
    
    bim_chrom =  genotypeData["chrom"]
    bim_pos =  genotypeData["pos"]
    bim_A1 = genotypeData["A1"]
    bim_A2 = genotypeData["A2"]

    del genotypeData ; gc.collect() # dont need this


    # 1. standardise data, this may remove columns, so we do this before we load PRS file
    if args.qc == 1 :
        start = time.time()

        qc_data = genoQC_all(M_orig, rsIds = irsIds, replaceMissing = _replaceMissing, minObserved = _minObserved, minMAF = _minMAF, minVariance = _minVariance) 
        rsIds_qc = qc_data["rsIds"] # save away the surviving SNP list that we have used 
        indicesToRemove = qc_data["indicesToRemove"]
        indicesKept = qc_data["indicesKept"]
        irsIds = rsIds_qc.tolist()
            
        del qc_data; gc.collect() # overwrite
            
        qc_data = removeList(M_orig, indicesToRemove, None, bim_chrom, bim_pos, bim_A1, bim_A2)
        M_orig = qc_data["X"]
        bim_chrom =  qc_data["bim_chrom"]
        bim_pos =  qc_data["bim_pos"]
        bim_A1 = qc_data["bim_A1"]
        bim_A2 = qc_data["bim_A2"]

        del qc_data; gc.collect() # overwrite
        end = time.time(); printElapsedTime(start,end, "QC took: ")
    else :
        M_orig[M_orig==-1] = 0 # have to make sure that the missing genotype is NOT encoded as -1, even when we don't perform QC
        print("Skipping internal QC", flush=True)

       
    # create a chr_pos lookup dictionary to be able to tell if the SNP we will load for the PRS is present at all
    bim_SNPs = list()
    bim_SNPs_dict = {}
    for i in range(len(bim_pos) ):   
        lookup = tuple(np.array([bim_chrom[i], bim_pos[i]], np.int32))
        bim_SNPs_dict[lookup] = [bim_A1[i], bim_A2[i] ] # save one lookup dict with the values of A1/A2 (to be used for matching SNPs)
        bim_SNPs.append(lookup )

    # load PRS file, if any
    if args.prs is not None : 
        PRS_SNPs_dict, PRS_SNPs_orig, PRS_betas_orig  = loadPRSFile(args.prs, bim_SNPs_dict)
        if suppressPrint == False : print("Loaded PRS for ",len(PRS_SNPs_orig) , " SNPs", flush=True)

        # get indices of SNPs in the genotype file, which we have PRS data for ( we keep the original order of the PLINK file)
        bim_SNP_indices = matchIndices(bim_SNPs, PRS_SNPs_dict)

        # find the indices of the SNPs in the PRS, now mapped to the subset of SNPs actually existing in both the PLINK and PRS files
        PRS_SNP_indices = matchIndices(bim_SNPs, PRS_SNPs_dict, True)
        PRS_SNPs = np.array(PRS_SNPs_orig)[PRS_SNP_indices]
        PRS_betas = PRS_betas_orig[PRS_SNP_indices]
        del PRS_betas_orig; del PRS_SNPs_orig

        #region EXPLANATION FOR IN-PLACE NUMPY ARRAY REARRANGEMENT
    ########################
    # example with the idx being a subset, ut where arr_sub has some non existent elements too
    #arr = np.array(["SNP1", "SNP2", "SNP3", "SNP4", "SNP5"])
    #arr_sub = np.array(["SNP2", "SNP5", "SNP3", "SNP6"])
    #arr_sub_betas = np.array([0.3, 0.7, 0.1])

    
    #dict_indices = {} # find the indices of the small array
    #for i in range(len(arr_sub) ): dict_indices[arr_sub[i]] = i # save antoher lookup dict that stores the indices 
    #idx_large = matchIndices(arr, dict_indices)
    #idx_large # array([1, 2, 4])


    #idx_small = matchIndices(arr, dict_indices, True)
    #idx_small # array([0, 2, 1])

    #arr[idx_large]            # array(['SNP2', 'SNP3', 'SNP5'], dtype='<U4')
    #arr_sub[idx_small]        # array(['SNP2', 'SNP3', 'SNP5'], dtype='<U4')
    #arr_sub_betas[idx_small]  # array([0.3, 0.1, 0.7])
    #bim_SNPs_subset = np.array(bim_SNPs)[bim_SNP_indices]
    #len(PRS_SNPs)
    #len(bim_SNPs_subset)
    #i=0
    #counter = 0
    #for i in range(len(PRS_SNPs) ):
    #    equal =  bim_SNPs_subset [ i ] [0] == PRS_SNPs[i][0] and bim_SNPs_subset [ i ] [1] == PRS_SNPs[i][1]
    #    if(equal == False) :  
    #        counter+=1
    #    print(i, " bim_SNPs: " , bim_SNPs_subset [i ] , " / PRS_SNPs: " , PRS_SNPs[i], "are equal: ", equal) 
    #print("not equal: " , counter, " / : " , len(PRS_SNPs) )
    ##########################################
    #endregion
    else : 
        bim_SNP_indices = np.asarray( range(len(bim_SNPs)) )
        PRS_betas = None


    pheno_dict = {}
    for i in range(len(indisWithPhenos) ): pheno_dict[indisWithPhenos[i]] = i

    # load covars, if any
    if args.covars_IDs is not None and (args.covars_cont is not None or args.covars_factor is not None):
        covarIDs = loadTextList(args.covars_IDs) 
        if args.covars_cont is not None : 
            covars_cont_orig= loadNumericFromTextFile(args.covars_cont) 
            if suppressPrint == False : print("Loaded continuous covars for ",covars_cont_orig.shape[0]," people and " , covars_cont_orig.shape[1] , " vars", flush=True)
        if args.covars_factor is not None : 
            covars_factor_orig= loadNumericFromTextFile(args.covars_factor, True) 
            if suppressPrint == False : print("Loaded factor covars for ",covars_factor_orig.shape[0]," people and " , covars_factor_orig.shape[1] , " vars", flush=True)

            # convert factors into one-hot encoded dummy design matrix. 
            # ie if we had [1,3,1] (with factor levels of 1, 3 and 2), this would get onverted to:
            #   [ 1, 0,   , 0,0,1,    ,1,0]
            # want to do this BEFORE we subset the data, in case we would loose factor level when excluding an individual
            # PROBLEM: if we are doing an inference run, then the inference subset may NOT have all of the factor levels, so in those cases we load those from disk to ensure that the same factor levels werere present both during training/testing)
            
            if args.inference == 0 :     # if we are running training, we assume that all factor levels are present in the data
                factorLevels = np.zeros( (covars_factor_orig.shape[1]), dtype = np.int32 ) # init with the right size
                for i in range(covars_factor_orig.shape[1] ): factorLevels[i] = int( np.max(covars_factor_orig[:,i])  ) # the number of factors is the max value
                 
                # write vector to disk
                writeVectorToDisk(args.out+ ".factorlevels",factorLevels, dataType ="int32")
            else : # for inerence runs we load the factor levels from disk
                factorLevels = loadVectorFromDisk(args.out+ ".factorlevels", dataType ="int32") 


            #region algo outline
            # eg: 3 factors, with max levels: [ 2, 3, 2]
            # [ 0,0,     0,0,0,    0,0]
            # we have indi with values: (1,3,2), 
            # should be coded as:
            # [ 1,0,     0,0,1,    0,1]
            # algo values:
            # i = 0: sumColsSoFar = 0, maxLevel = 2
            # [j, 0 + 1-1 ] = 1  //  [ 1,0,
            # i = 1: sumColsSoFar = 2, maxLevel = 3
            # [j, 2 + 3-1 ] = 1  //  [ 1,0,     0,0,1,
            # i = 2: sumColsSoFar = 5, maxLevel = 2
            # [j, 5 + 2-1 ] = 1  //  [ 1,0,     0,0,1,    0,1]
            #endregion
            # transcode the factor levsl into the dummy design matrix  
            covars_factor_orig_dummy = np.zeros( (covars_factor_orig.shape[0], sum(factorLevels)), dtype = np.int32 ) # init it with the right size
            sumColsSoFar = 0
            #i= 1; j = 0
            for i in range( len(factorLevels) ) : # go through each of the original columns/factor levels
                maxLevel = factorLevels[i]
                for j in range( covars_factor_orig.shape[0] ) : # go through each indi
                    covars_factor_orig_dummy[j,  sumColsSoFar+ int( covars_factor_orig[j,i] )-1  ] = 1 # 1 hot, with the index coming from the actual number, offset by the sumColsSoFar

                sumColsSoFar = sumColsSoFar + maxLevel
            if suppressPrint == False : print("transcoded factor covars as one hot for  shape",covars_factor_orig_dummy.shape[0]," people and " , covars_factor_orig_dummy.shape[1], flush=True)
            covars_factor_orig = covars_factor_orig_dummy
            del covars_factor_orig_dummy
        
        # get keeplist of individuals, which is the intersection of the phenos, covars and the actual genotype files
        covar_dict = {}
        for i in range(len(covarIDs) ): covar_dict[covarIDs[i]] = i

        Indis_in_plink_indices = matchIndices2(IDs[0], covar_dict, pheno_dict) # will need to use the 2 dict version, as we only want to keep indis, who have both pheno and covar files
    else :
        Indis_in_plink_indices = matchIndices(IDs[0], pheno_dict)


    Indis_in_plink = IDs[0][Indis_in_plink_indices] # get the final list of indis, in the order that they will be in the subset genotype file
 
    # find the indices of the extant indis in the covars/phenos
    # subset each (the keeplist indices kept must match each modality)
    if args.covars_IDs is not None and (args.covars_cont is not None or args.covars_factor is not None):
        covars_indices_in_bim = matchIndices(Indis_in_plink, covar_dict, True)
        #covarIDs_keeplist = np.array(covarIDs)[ covars_indices_in_bim ]
        if args.covars_cont is not None : 
            covars_cont = covars_cont_orig[covars_indices_in_bim,:]
            del covars_cont_orig
        else : covars_cont = None

        if args.covars_factor is not None : 
            covars_factor = covars_factor_orig[covars_indices_in_bim,:]
            del covars_factor_orig
        else : covars_factor = None
    else :
        covars_all = None
        covars_cont = None
        covars_factor = None

    # subset phenos to the subset of indis that had all covariates
    phenos_indices_in_bim = matchIndices(Indis_in_plink, pheno_dict, True)
    #phenos_keeplist = np.array(indisWithPhenos)[ phenos_indices_in_bim ]

    y_keep = y[ phenos_indices_in_bim ]
    y = y_keep
    del y_keep


    # subset the genotype matrix
    #M_keep = M[phenos_indices_in_bim,bim_SNP_indices] # cannot simply use 2 indices, as I will get error: IndexError: shape mismatch: indexing arrays could not be broadcast together with shapes
    M = M_orig[np.ix_(phenos_indices_in_bim,bim_SNP_indices)] # have to use np.ix_ : https://stackoverflow.com/questions/35607818/index-a-2d-numpy-array-with-2-lists-of-indices
    del M_orig

    np.random.seed(args.randomSeed)
    random.seed(args.randomSeed)  
    if args.inference == 0 :     # if we are running training

        # load the validation split file, we need this, as we only want to standardize over the training set: and not the training+ validation together: https://blog.slavv.com/37-reasons-why-your-neural-network-is-not-working-4020854bd607
        validIDs = loadTextList(args.validSet) # if we are in analysis mode, then a validation set MUST have been provided so it is safe to call this
        
        # create subset of indices for the training and validation 
        training_indices, valid_indices = validTrainSplit(validIDs, Indis_in_plink)

        len_M_validation = len(valid_indices)
        print( str( len(training_indices) ) + " training set and " + str( len(valid_indices) )+ " valid set", flush = True)
        #numIndividuals_valid = len(valid_indices)

        IDs_validation =  Indis_in_plink[valid_indices]

        start = time.time()
        
        if _standardiseData :
            # if  Training (IE not inference), we standardise Y right on the spot from the training data
            y_training = y[training_indices]
            y_training, y_mean, y_sd = zscore(y_training) # zscore it so that Beta -> h2 computations work    

            y_mean = y_mean[0] ; y_sd = y_sd[0]
            writeYMeanSds(args.out, y_mean,y_sd ) # write the y_mean/valid to disk so that it could be loaded for test sets
            y = y_standardise(y,y_mean,y_sd) # standardise the whole y vector but from the mean/sd,  from the training y only
            # do NOT standardise genotypes. Only continuos covars
            del y_training;
            #M, mns, sstd = standardise_Genotypes(M) ; gc.collect()
            if args.covars_cont is not None :
                covars_cont_training = covars_cont[training_indices] # again, use just the training set, not the whole set to compute the mean/sd
                covars_cont_training, mns, sstd = zscore(covars_cont_training) ; 
                covars_cont -= mns
                covars_cont /= sstd

                #write training data means / stds to disk so that we could use those for inference runs later
                print("writing means/stds to disk with datatype: "  + str(sstd.dtype))
                print("sstd shape is: " + str(sstd.shape) + " / mns shape: " + str(mns.shape))
                knet_IO.writeVectorToDisk( args.out + "data_mns" , mns, mns.dtype)  
                knet_IO.writeVectorToDisk( args.out + "data_sstd" , sstd, sstd.dtype)  

                del covars_cont_training; gc.collect()

        end = time.time(); printElapsedTime(start,end, "standardising data took: ")
        print("After standardising, training data in MBs is: ",getSizeInMBs(M) )

        
        # apply oversampling logic for case control phenos, if requested (this does not cause any increase in memory usage) (do this AFTER we have standardised, as standardisation depends on SD, which could be changed if we oversampled)
        if args.oversampling is not None and cc is True :     
            caseIndices = list()
            controlIndices = list()

            # find the cases/controls indices, it is assumed the cases are coded as '1'
            for i in range(len(training_indices) ): # want to loop over the training indices subset not directly on the whole of y!!
                index = training_indices[i]
                if y[index][1] == 1 : caseIndices.append(index)  # y[i][ 1 ], as the cases are the 2nd column in a one-hot data shape
                else : controlIndices.append(index)

            print("num cases : " + str(len(caseIndices)) + " and num controls: " + str(len(controlIndices)), flush=True)
            if len(caseIndices) != len(controlIndices) :
                print("applying oversampling logic as these were not equal...", flush=True)
                
                # sample with replacement so that the shorter list will be as long as the longer
                if len(controlIndices) < len(caseIndices) : controlIndices = choices(controlIndices, k=len(caseIndices))  
                else : caseIndices = choices(caseIndices, k=len(controlIndices))

                # to create a balanced training set, we just combine the two lists, and overwrite the training_indices with it. This will have some duplicate elements, but that is OK, we are not duplicating the memory footprint of the plink genotype data as we are only duplicating the lookups
                training_indices = [*controlIndices, *caseIndices] # concat the 2 lists, and convert it to a numpy array
                training_indices = np.array(training_indices)

        # Shuffle data before producing the minibatches to avoid having all-case or all-control minibatches
        # only shuffle order for analysis runs, as for inference we want to keep the original
        start = time.time()
  
        # we only shuffle the training indices, not the actual data, as that is never directly used, we always go through the training indices to extract a minibatch etc
        indices = np.asarray( range(len(training_indices)) ) # is used for storting
        random.shuffle(indices)
        training_indices = training_indices[indices]

        end = time.time(); printElapsedTime(start,end, "shuffling data took: ")

    else : # inference run
        print("Inference data QC", flush=True)
        if args.snpIndices is not None :
            indicesToKeep = knet_IO.loadIndices(args.snpIndices)
            M = M[:,indicesToKeep]
        len_M_validation = 0 
        start = time.time()
        if _standardiseData :
            # if  inference, we standardise Y from what is on the disk
            y_mean, y_sd = loadYMeanSds(args.out+ ".ymeansd")
            
            if y_mean is None : # if we were unable to load the mean/sd from disk (like if we did not calculate it)
                print("y mean/sd was not pre calculated previously, we do it on the fly")
                y, y_mean, y_sd = zscore(y) # zscore it so that Beta -> h2 computations work    
                y_mean = y_mean[0] ; y_sd = y_sd[0]
            else :
                y = y_standardise(y,y_mean,y_sd)
            
        M[M==-1] = 0  # have to make sure that the missing genotype is NOT encoded as -1, even when we don't perform QC

        if _standardiseData and args.covars_cont is not None  : # standardise the continous covars if there were any (NOT THE GENOTYPES!)
            # standardise covar data too
            mns  = loadVectorFromDisk( args.mns  , 'float32')  # these are always float32 even in 64 runs
            sstd = loadVectorFromDisk( args.sstd , 'float32')  
            covars_cont -= mns
            covars_cont /= sstd

        training_indices =  np.asarray( range(len(M)) ) # training indices is simply the unshuffled 1:numIndis
        end = time.time(); printElapsedTime(start,end, "standardising data via loaded params took: ")

        
    # merge continuous and factor covars (after the former has been standardised)
   
    if covars_cont is not None or covars_factor is not None  :
        if covars_cont is None : covars_all = covars_factor # if we only have one or the other we just use that as the 'all'
        elif covars_factor is None : covars_all = covars_cont
        else : # if we have both we have to merge them
            covars_all = np.zeros( (covars_factor.shape[0],covars_factor.shape[1] + covars_cont.shape[1] ), dtype = np.float32 ) # create blank 2D array that could accommodate both covars data
            covars_all[:,0:covars_factor.shape[1]] = covars_factor # paste in the factors data in the first part
            covars_all[:,covars_factor.shape[1]:covars_factor.shape[1] + covars_cont.shape[1]] = covars_cont # paste in the continous data in the second half

    del covars_cont; del covars_factor; gc.collect(); # these were inited to be None, so it is safe to call delete


    # 2. create minibatch list (after we have subset and split training/validation)
    #numIndividuals = M.shape[0] 
    global numSNPs; numSNPs = M.shape[1] # need this global, as in the convolution scenario, batch.shape[1] will always be 1
    num_y_classes = 1 # how many columns are there, IE how many classes 
    len_M = len(training_indices)
    global numCovars; 
    numCovars = 0
    if covars_all is not None  : numCovars= covars_all.shape[1]


    # 3. create minibatches (cant do this in a function, as that would make it unavoidable to use 2x the RAM)
    training_indices_orig = np.copy(training_indices) # create backup of original, as we will eliminate it via minibatch creation
    valid_indices_orig = None
    startTime = time.time()
    start = 0
    minibatch_size =  args.batch_size #M.shape[0]  # 64   #minibatch_size = 128
    if args.batch_size == 0 : minibatch_size = len(M)
    num_batches = len_M // minibatch_size
    num_batches = max(1,num_batches) # make sure to have at least 1 minibatch...
    end = minibatch_size
    train_minibatchIndices = list()

    for i in range(num_batches) :  # # do this in a more RAM efficient way: keep deleting the bits from the original matrix to free up space as we go along otherwise this step would double the RAM requirements temporarily
        train_minibatchIndices.append(training_indices[0:minibatch_size]  )
        training_indices = training_indices[minibatch_size:len(training_indices)]

        print("adding batch " + str(i)  + ", minibatch size: " + str(minibatch_size) + " / num left in pool: " + str(len(training_indices))  )
        gc.collect()

    remainderBatch = training_indices


    remainderBatch_valid = None   
    if args.validSet is not None : # if there is a separate validation set (for inference runs we won't have this as there is only a test set (which is confusingly stored as the training set)
        len_M_validation = len(valid_indices) 
        valid_indices_orig = np.copy(valid_indices) 
        if args.batch_size == 0 : minibatch_size = len(valid_indices)
        num_batches = len(valid_indices) // minibatch_size
        num_batches = max(1,num_batches) # make sure to have at least 1 minibatch....
        test_minibatchIndices = list()
        print("len_M_validation is: " + str(len_M_validation)  + ", minibatch size: " + str(minibatch_size)  + " args.batch_size: " + str(args.batch_size) + " num_batches is: " + str(num_batches))
        start = 0
        end = minibatch_size
        for i in range(num_batches) :
            #test_GWAS.append(M_validation[start:end]  )
            test_minibatchIndices.append(valid_indices[0:minibatch_size]  )
            valid_indices = valid_indices[minibatch_size:len(valid_indices)]
            print("adding batch " + str(i)  + " , start/end: " + str(start) + "/" + str(end)  )
            start = end
            end += minibatch_size  
        remainderBatch_valid = valid_indices
        #print("First minibatch is: ", test_GWAS[0][0,])
        # del M_validation; gc.collect() # free up memory, cant do this as we need this for the PRS calculation....

    else :
        test_minibatchIndices = None 
    end = time.time(); printElapsedTime(startTime,end, "creating minibatches took: ")    
    #print("Last minibatch is: ", train_GWAS[-1][0,])
    #print("Last minibatch As Valid is: ", test_GWAS[-1][0,])   

    # scale the delta by minibatch_size, if we dont have minibatches
    #ratio = float(minibatch_size) / numIndividuals # this is 1 if there are no minibatches
    #print("orig L2 Regularizer : " + str(args.hidl2) + " minibatches scaled to " + str(hiddenShrinkage * ratio) )
    #hiddenShrinkage *= ratio


    #endregion
    
    ###########################################################################
    # II) Model Build: find or load the best hyper param settings and build a model
    ##########
    # setup device
    global device
    device = torch.device("cuda:0" if torch.cuda.is_available() and args.gpu ==1 else "cpu")
    args.accumulation_steps = max(1,args.gradient_batch_size // args.batch_size ) # number of accumulation steps is based on the effective batch size we want for the gradient calculation ( do this before upscalign batch size)
    print("gradient batch_size:", args.gradient_batch_size, " / batch_size:", args.batch_size, "/ accumulation_steps worked out to be", args.accumulation_steps)


    # get untrained random baseline accuracy
    train_chanceLevel_loss = meanLevelLoss(y, training_indices_orig) ; 
    if valid_indices_orig is not None: valid_chanceLevel_loss = meanLevelLoss(y, valid_indices_orig)
    else : valid_chanceLevel_loss = -1
    line="Trainset len/batches: " +str(len(y) ) + "/"+ str(len(train_minibatchIndices) )  + " | chance level prediction training Loss: " + str(round( train_chanceLevel_loss,3))   + " | validation loss: " + str(round(valid_chanceLevel_loss ,3) )   ; print( line)
    with open(args.out+"_chance", "w") as file: file.write(line  + "\n")

    # 4a: if hyoperopt was enabled then determine the best hyperparam settings automatically
    if args.hyperopt != 0 and args.inference == 0 : # can only happen if it is not inference run
        startups = max(5,int(args.hyperopt / 4) )
        print("determining best params via hyperopt for num trials: ", args.hyperopt, " / startup jobs: ", startups)
        start = time.time()
        best_pars = optimize_model_pytorch(device, args, y, M, train_minibatchIndices, PRS_betas, covars_all, test_minibatchIndices, out_folder = args.out +"hyperopt/", startupJobs = startups, maxevals = args.hyperopt)     #5, maxevals = 20)   
        writeKNeT_bestPars(args.out ,best_pars)    
        #best_pars = loadKNeT_bestPars(args.out)
        args = mergeArgsAndParams(args,best_pars)
        torch.cuda.empty_cache()
        end = time.time(); printElapsedTime(startTime,end, "hyperopt took: ")

    ###########################################################################
    # III) Train model:
    ###################
    # 4b. build final model 
    # if a linear activation function was selected, then we enable batchnorm. Otherwise this would downward bias the non-linear estimates,
    # as there the 'non-linear' model would still be linear, but without batchnorm which would destroy the non-linear accuracies
    if args.hidAct == 0 : 
        args.bnorm = 1
        print("as a linear activation was selected, we enable bnorm even for the nonlinear model")
            
    model = build_model(args, numSNPs,  numCovars, num_y_classes) # Don't pre-allocate memory; allocate as-needed ??
    model_diagnostics(model,args)
    model = modelToDevice(model,args,device) # move model to device (CPU or GPU)

    # 6a. Analysis: train model
    if args.inference == 0 :
        print("Training  moden starts", flush = True)
        postFix = ""
        if modelAlreadyRun (args, postFix = postFix) == False :
            trainModel(args, model, train_minibatchIndices,remainderBatch_valid, test_minibatchIndices, M, IDs_validation, PRS_betas,  covars_all, len_M_validation, postFix = postFix, y_mean = y_mean, y_sd = y_sd, indicesKept = indicesKept, irsIds = irsIds)  

        postFix = linearPostFix
        if modelAlreadyRun (args, postFix = postFix) == False :
        # for the linear model we need to enable bnorm and retrain from scratch, as switching off the SELU would remove the standardisation feature as well, which would make linear models unfairly poor
            print("RETRAINING A LINEAR MODEL FROM SCRATCH")
            args.bnorm = 1 # bnorm is always enabled for linear models
            args.hidAct = 0;  # switch activation off
            model = build_model(args, numSNPs,  numCovars, num_y_classes) # Don't pre-allocate memory; allocate as-needed ??
            model = modelToDevice(model,args,device) # move model to device (CPU or GPU)
            trainModel(args, model, train_minibatchIndices,remainderBatch_valid, test_minibatchIndices, M, IDs_validation, PRS_betas,  covars_all, len_M_validation, postFix = postFix, y_mean = y_mean, y_sd = y_sd, indicesKept = None, irsIds = None)  

    # 6. b analysis: inference we build polygenic risk scores
    else :
        print("Inference Run", flush = True)
        if torch.cuda.is_available() : model.load_state_dict(torch.load(args.loadWeights))
        else :  model.load_state_dict(torch.load(args.loadWeights,  map_location='cpu'))

        profileName= "yhat_TEST.txt"
        rSQoutName ="KNET_PRS_TEST"
            
        # the Train set here will refer to the TEST set
        producePRS(model,remainderBatch, train_minibatchIndices, M, IDs    , PRS_betas,  covars_all, len_M , args.out + profileName, args.out + "FIDs_TEST.txt", y, args.out + rSQoutName, y_mean, y_sd)  # (model,args,remainderBatch, train_minibatchIndices, IDs, len_M , args.out + profileName, args.out + "FIDs_TEST.txt", y, args.out + rSQoutName, y_mean, y_sd)
    
        # produce _noAct (this is the same as the non-linear model, we just switch off the activation layers)
        model = getInferenceModelWithoutActivation(model)  #if args.linearInference == 1 :
        profileName = "yhat_TEST_noAct.txt" 
        rSQoutName ="KNET_PRS_TEST_noAct"
        producePRS(model,remainderBatch, train_minibatchIndices, M, IDs    , PRS_betas,  covars_all, len_M , args.out + profileName, args.out + "FIDs_TEST.txt", y, args.out + rSQoutName, y_mean, y_sd)
       
          
        # Produce Retrained _noAct
        linearModelWeights = args.loadWeights+ linearPostFix
        my_file = Path( linearModelWeights)
        if my_file.is_file():
            args.bnorm = 1 # bnorm is always enabled for linear models
            args.hidAct = 0;  # switch activation off
            print("a retrained linear model also exist, we produce PRS for that too")
            model = build_model(args, numSNPs,  numCovars, num_y_classes) # Don't pre-allocate memory; allocate as-needed ??
            model = modelToDevice(model,args,device) # move model to device (CPU or GPU)

            if torch.cuda.is_available() : model.load_state_dict(torch.load(linearModelWeights))
            else :  model.load_state_dict(torch.load(linearModelWeights,  map_location='cpu'))
            profileName = "yhat_TEST_noAct_retrain.txt" 
            rSQoutName ="KNET_PRS_TEST_noAct_retrain"

            producePRS(model,remainderBatch, train_minibatchIndices, M, IDs    , PRS_betas,  covars_all, len_M , args.out + profileName, args.out + "FIDs_TEST.txt", y, args.out + rSQoutName, y_mean, y_sd)
        else: 
            print("retrained linear model does not exist")



# trains and evaluates a model
def trainModel(args, model, train_minibatchIndices,remainderBatch_valid, test_minibatchIndices, M, IDs_validation, PRS_betas,  covars_all, len_M_validation, postFix = "", y_mean = None, y_sd = None, indicesKept = None, irsIds = None) : 
    saveModelLocation = None
    if args.saveWeights is not None : saveModelLocation = args.saveWeights + postFix
    # attempt to resume previous model, if any
    t, training_hasntmprovedNumIts, results, optimizerState =  attempResumeModel(args, model,device, saveModelLocation)

    start = time.time()                                                                                                                                                            
    results = learn(model, device, args, y, M, train_minibatchIndices, PRS_betas=PRS_betas , covars_all =covars_all, test_minibatchIndices=test_minibatchIndices, saveModelLocation = saveModelLocation, epochMaxImproveThreshold = args.epochMaxImproveThreshold, learnRate =args.learnRate, momentum = args.momentum, half = args.half, decayRate = args.LRdecay, accumulation_steps =args.accumulation_steps, l2= args.l2, optimizerChoice = args.optimizer, pheno_is_binary = pheno_is_binary, t = t, training_hasntmprovedNumIts = training_hasntmprovedNumIts, results = results, optimizerState = optimizerState) # (model, device, args, y, M, train_minibatchIndices, PRS_betas=PRS_betas , covars_all =covars_all, test_minibatchIndices=test_minibatchIndices, eval_train=True, eval_test=True, eval_freq = args.evalFreq, decayEnabled = False)    
    end = time.time(); printElapsedTime(start,end, "training model took: ")
    results_its = results["results"]  
      
    # 3)   write out the best epoch (this is 0 based)
    with open( args.out + postFix + best_epoch_string, "w") as file: file.write("best_epoch=" + str(results['results']['lowestLoss_epoch']) ) # write out the early stop epoch used


    # 4) Save model params
    if args.saveWeights is not None : 
        lowestLoss_epoch =results['results']['lowestLoss_epoch']
        for i in range( len(results['results']["epochs"]) ) :
            epoch = results['results']["epochs"][i]
            
            if epoch == lowestLoss_epoch : # rename the best model 
                print("highest performing model at epoch",lowestLoss_epoch, " saved to:",args.saveWeights+ postFix)
                os.rename(args.saveWeights+ postFix + str(lowestLoss_epoch), args.saveWeights+ postFix)
                os.rename(args.saveWeights+ postFix + str(lowestLoss_epoch)+ _results , args.saveWeights + _results + postFix)
                os.rename(args.saveWeights+ postFix + str(lowestLoss_epoch)+ _optimizer, args.saveWeights+ _optimizer+ postFix)
                os.rename(args.saveWeights+ postFix + str(lowestLoss_epoch)+ _stats , args.saveWeights+ _stats + postFix)
            elif epoch > 0 : # delete the rest, except model0, the one with no training as a reference
                os.remove(args.saveWeights+ postFix + str(epoch))
                os.remove(args.saveWeights+ postFix + str(epoch)+ _results)
                os.remove(args.saveWeights+ postFix + str(epoch)+ _optimizer)
                os.remove(args.saveWeights+ postFix + str(epoch)+ _stats)


    ###########################################################################
    # IV) Output diagnostics:
    #########################
    print("(IV) Output", flush = True)

    # 1) write out summary results of the training iterations
    fileName = args.out + postFix + "nn_results.txt"
    with open(fileName, "w") as file:      
        line = "epochs"
        if "train_loss" in results_its: line = line + "\t" + "train_loss"
        if "valid_loss" in results_its: line = line + "\t" + "valid_loss"

        file.write(line  + "\n")
         
        for i in range( len(results_its["epochs"])  ):
            line = str(results_its["epochs"][i]) 
            if "train_loss" in results_its: line = line + "\t" + str(results_its["train_loss"][i])
            if "valid_loss" in results_its: line = line + "\t" + str(results_its["valid_loss"][i])

            file.write(line + "\n")   


    # 2) generate learning curve plots of the results
    if len(results_its["epochs"]) > 0 :
        metric = "MSE"
        if pheno_is_binary : metric ="BCE"
        plotNNtraining(results_its, args.out + "trainplot_loss"+postFix, training_metric="train_loss", valid_metric="valid_loss", acc_measure = metric)

    # 3a) output validation accuracy: of the FINAL model
    outputAccuracy(args, model,remainderBatch_valid, test_minibatchIndices, M, IDs_validation, PRS_betas,  covars_all, len_M_validation, postFix = "_final", y_mean = y_mean, y_sd = y_sd) 


    # 3b) output the validation accuracy for the BEST model
    modelWeightsToload = args.saveWeights+ postFix
    model = build_model(args, device, mapModelToDevice = False) # do NOT parallelise it yet or move it to the GPU
    model.load_state_dict(torch.load(modelWeightsToload,  map_location='cpu')) # load it to the CPU first
    model = modelToDevice(model,args,device)
    outputAccuracy(args, model,remainderBatch_valid, test_minibatchIndices, M, IDs_validation, PRS_betas,  covars_all, len_M_validation, postFix = "_best", y_mean = y_mean, y_sd = y_sd) 
    
    # write out the SNPs that were used for the analysis
    if irsIds is not None :
        fileName = args.out + "nn_SNPs.txt"
        with open(fileName, "w") as file: 
            for i in range( len(irsIds)  ):
                file.write(irsIds[i]  + "\n")
        
    # write out the indices of the original dataset's coordinates for convenience
    if indicesKept is not None: # in case we skipped QC
        fileName = args.out + "nn_SNPs_indices.txt"
        with open(fileName, "w") as file: 
            for i in range( len(indicesKept)  ):
                file.write( str(indicesKept[i])  + "\n")    
             

#endregion        
         
    
###############################################################################
# Helper functions
#region Helper functions
###############################################################################

# tries to load and resume an existing model
def attempResumeModel(args, model,device, saveModelLocation) :
    # 2) attempt to load model weights if a previous run hasn't finished
    
    # keep looping to check if a model for a certain epoch exists:
    lastExistingEpoch = -1
    t = 0
    existingModel= saveModelLocation+ str(t) + _results
    
    while os.path.exists(existingModel ) :
        lastExistingEpoch = t
        t = t +1
        existingModel= saveModelLocation+ str(t) + _results
        
    if lastExistingEpoch != -1 :
        print("resuming previous model from epoch: ", lastExistingEpoch)
        oldModelLoc = saveModelLocation + str(lastExistingEpoch)
        
        # load the stats
        results, t, training_hasntmprovedNumIts = loadEpochStats(oldModelLoc)
        
        # load the model
        if torch.cuda.is_available() and args.gpu > 0 : 
            
            if args.device != -1 : # if a specific device Id was requested, we try to map weights there
                model.load_state_dict(torch.load(oldModelLoc,  map_location="cuda:"+str(args.device)))
            else :model.load_state_dict(torch.load(oldModelLoc)) # otherwise attempt reconstruct the model to same devices as they were trained on, IE for multi GPUs
            
        else :  model.load_state_dict(torch.load(oldModelLoc,  map_location='cpu'))
        
        optimizerState = torch.load(oldModelLoc+ _optimizer)
            
        return (t+1), training_hasntmprovedNumIts, results, optimizerState # +1 as we ALREADY have t, so we want the next one
        
    # if there wasn't any previous model ,we return the defaults
    else : return 0,0,None  , None



# generates yhat output of a model, and if it was binary pheno, then also ROC/AUCs too
def outputAccuracy(args, model,remainderBatch_, minibatchIndices, M, IDs_, PRS_betas,  covars_all, len_M_, postFix ="", y_mean = None, y_sd = None) :
    # predict validation accuracy 
    yhat, y = producePRS(model,remainderBatch_, minibatchIndices, M, IDs_, PRS_betas,  covars_all, len_M_ , args.out + postFix + "yhat.txt", args.out  + postFix+ "FIDs.txt", y, args.out  + postFix+ "KNET_PRS", y_mean, y_sd)   # (model,args,remainderBatch_valid, test_minibatchIndices, IDs_validation, len_M_validation , args.out + "yhat.txt", args.out + "FIDs.txt", y_validation, args.out + "KNET_PRS", y_mean, y_sd)         
            
    # for binary classifications we also want to produce a AUC ROC
    if pheno_is_binary :
        # as ROC
        false_positive_rate, true_positive_rate, thresholds = roc_curve(y[:,0], yhat[:,0])
        roc_auc = auc(false_positive_rate, true_positive_rate)
        
        plotROC_or_PrecisionRecall(args.out + postFix + "_ROC",false_positive_rate, true_positive_rate, roc_auc)
        
        # as Precision-recall
        # no Skill: a model that would predict the majority class 
        counts = np.bincount(y[:,0].astype(dtype=int)) # find most common element: https://stackoverflow.com/questions/6252280/find-the-most-frequent-number-in-a-numpy-array
        majorityClass = np.argmax(counts)
        yhat_noskill = np.full(y.shape, majorityClass) # create a yhat prediction of the majority class for each sample
        precision, recall, thresholds = precision_recall_curve(y[:,0], yhat_noskill[:,0])
        no_skill = auc(recall, precision)

        precision, recall, thresholds = precision_recall_curve(y[:,0], yhat[:,0])
        auc_precisionRecall = auc(recall, precision)

        plotROC_or_PrecisionRecall(args.out+ postFix + "_PR",recall, precision, auc_precisionRecall, False, no_skill)
       
        with open(args.out + postFix + "model_auc.acc", "w") as file:    
            file.write("ROC_AUC:\t" + str(roc_auc) + "\n" ) 
            file.write("PR_AUC:\t" + str(auc_precisionRecall) + "\n" ) 
    

# moves the model to the appropriate device (CPU or a GPU or even multiple GPUs), depending on what was requested
def modelToDevice(model,args,device, suppressPrint = False) :
    # data parallelism: https://pytorch.org/tutorials/beginner/blitz/data_parallel_tutorial.html
    if torch.cuda.device_count() > 1 and args.gpu  > 1  : # do NOT use dataparallel GPU for inference runs as the hooks don't work: # these may not work on DataParallel models ??: https://pytorch.org/docs/0.3.1/nn.html?highlight=hook#dataparallel-layers-multi-gpu-distributed
        if suppressPrint == False : print("Multiple GPUs requested", args.gpu)
        gpuIDs = np.array( range(args.gpu) ).tolist() # determine the number of GPUs to use from the number supplied
        #gpuIDs = [2,3,4,5,6,7] # the first of these MUST be the 'device', IE the GPU that stores the master copy otherwise pytorch will error out
        model = nn.DataParallel(model , device_ids=gpuIDs)
        
    if suppressPrint == False : print(" args.gpu:",  args.gpu, "torch.cuda.device_count():", torch.cuda.device_count())
    model.to(device)
    return(model)


def loadStringList(outFile) :  # used by the oversampling logic to load in a list of cases
    items = list()
    with open(outFile, "r") as file: 
        for i in file:
            itmp = i.rstrip().split()
            items.append(itmp)
            
    items = np.array(items)

    return( items)   
    

def getInferenceModelWithoutActivation(model) : # produces an identical model, but without the activation layers (IE to get a linear predictor)
    print("switching off activations for linear infernece")
    model = getModel(model)
    origLayers = list(model)
    subsetModelLayers = list()
    for i in range(len(origLayers)) :
        if isLayerActivation(origLayers[i]) == False or i == ( len(origLayers) -1 ) : subsetModelLayers.append(origLayers[i]) # we add the last layer, even if that is activation, as that is needed to get the right shaped output
    modelInference = nn.Sequential(*subsetModelLayers) # create a subset model of onl y
    #modelInference.eval() # otherwise dropout and other layers would act up

    return(modelInference)


# miniBatches = train_minibatchIndices
def producePRS(model,remainderBatch      , miniBatches          , M, IndiIDs       , PRS_betas,  covars_all, len_total       , outLoc_yhat, outLoc_FIDs, ytrue, outLoc_PRS, y_mean = None, y_sd = None) :
    model_training_orig = model.training # otherwise dropout and other layers would act up
    setModelMode(model, False)
    global device
    # write final predictions out
    yhats = list()
    totalSofar= 0
    ytrue_subset = list() # need a separate array, as when we compute the correlation we want it to be for the subset we have actually produced results for, rather than the full y, which potentially includes both training/validation
    with torch.no_grad(): # reduce memory usage and speed up computations but you won’t be able to backprop; source: https://discuss.pytorch.org/t/model-eval-vs-with-torch-no-grad/19615
        for i in range(len(miniBatches)) : # loop through all minbatches
            totalSofar += len(miniBatches[i])
            b_labels, b_data = getMiniBatchData(ytrue, M, miniBatches[i], PRS_betas,  covars_all)
            b_data = torch.from_numpy(b_labels).to(device)  # b_data = torch.from_numpy(miniBatches[i]).to(device)
            yhats.append( model(b_data).detach().cpu().numpy() )
            ytrue_subset.append( ytrue[ miniBatches[i]] )
       
        if totalSofar < len_total :
            print("minibatches did not cover all training samples, so we create last batch out of the remainders")
            #lastBatch_X = origData[totalSofar:len_total]
            b_labels, b_data = getMiniBatchData(ytrue, M, remainderBatch, PRS_betas,  covars_all)
            b_data = torch.from_numpy(b_data).to(device)  # b_data = torch.from_numpy(remainderBatch).to(device)
            yhats.append( model(b_data).detach().cpu().numpy() )
            ytrue_subset.append( ytrue[remainderBatch] )
    
    yhat_all = np.concatenate(yhats)
    ytrue_subset = np.concatenate(ytrue_subset)
    print("after merging, we have yhat predictions for : " + str(len(yhat_all)) + " samples", flush=True)
    setModelMode(model, model_training_orig)  # reset model into traning mode
   
    # write out the final r^2
    yhat_all += EPSILON # for numerical stability




    rSQ = np.corrcoef( ytrue_subset, yhat_all, rowvar=0)[1,0]**2      
       
    with open(outLoc_PRS, "w") as file: 
            file.write(str(rSQ) ) 

    # transform both y and y hat back onto the raw scale, if this was requested
    if y_mean is not None :
        yhat_all = y_origscale(yhat_all, y_mean, y_sd)
        ytrue_subset = y_origscale(ytrue_subset,y_mean, y_sd) 
    
    fileName = outLoc_yhat
    with open(fileName, "w") as file:
        file.write("Profile"  + "\n")
        
        for i in range(yhat_all.shape[0]) :
            line = str(yhat_all[i][0] )
            for j in range(1, len(yhat_all[i]) ):
                line = line + "\t" + str(yhat_all[i][j] )              
            file.write( line +  "\n")   # file.write(  ( str(yhat[i])[2:-1] ).replace("  ", " ").replace(" ", "\t") +  "\n")

    # also write out the FID / IIDs in the same order, just as a sanity check (compare this against the .fam files)
    fileName = outLoc_FIDs
    with open(fileName, "w") as file:
        file.write("FID" + "\t" + "IID" + "\n")

        for i in range( len(IndiIDs[0]) ) :
            line = IndiIDs[0][i] + "\t" + IndiIDs[1][i]
            file.write( line +  "\n") 

    return yhat_all, ytrue_subset
#endregion     

###############################################################################
# Inference plotting
#region Inference plotting
###############################################################################
targetLayerActivations = None

def get_targetLayerActivations():
    
    def hook(model, input, output):
        global targetLayerActivations
        #if model.training : # ONLY add this if the model is in training mode, bad idea, as we only 
        targetLayerActivations = output  # output.detach()

    return hook

def produceActivation (model, device,targetLayerIndex,artificial_SNP_data):

    artificial_SNP_data = torch.from_numpy(artificial_SNP_data).to(device)
     
    origLayers = list(model)
    # targetLayerIndex = findPrevNextActivationLayer(model,targetLayerIndex, startIndexOK = True) # dont do this, we should assume that the correct layer was already chosen outside
    if targetLayerIndex < 0 : targetLayerIndex = len(model) +targetLayerIndex
    lastLayerSlice = targetLayerIndex +1 # , +1 as this is a slice, that is EXCLUSIVE, IE 0:len
    subsetModelLayers = origLayers[0:lastLayerSlice]
    modelInference = nn.Sequential(*subsetModelLayers) # create a subset model of onl y
    modelInference.eval() # otherwise dropout and other layers would act up
    
    # setup hook for FP, to capture intermediate activation
    activationHook = modelInference[targetLayerIndex].register_forward_hook(get_targetLayerActivations()) 
    
    modelInference(artificial_SNP_data) # FP the data
    
    global targetLayerActivations
 
    targetLayerActivations = targetLayerActivations.detach().cpu().numpy() # obtain the activation, this would have been set via the hook

    # remove hook
    activationHook.remove()         
    return(targetLayerActivations)
    
    
# targetLayerIndex has to be CORRECT
def produceActivationPlot(model,device,interactionsToTest,numNeurons, targetLayerIndex, totalNumSNPs, plotFolder_used, strength = 1, outFileStem = "true", normalized = None , scale = [8,8], rSQ = -1, doAll = False, subtractNullActivations = True ) :
    offset = 1
    if doAll == False : offset = 0
    activationMap = np.zeros( (len(interactionsToTest) +offset,numNeurons) ) # the activation map will have 1 row for each SNP, (AND +1 for when all are active), and one column for each neuron
    SNPlabels = list()
    neuronlabels = list()
    neuronlabels = list(range(numNeurons)) 
    neuronlabels =[x+1 for x in neuronlabels]
    #nullActivations.shape
    nullActivations = produceActivation(model,device,targetLayerIndex, np.zeros( (1,totalNumSNPs) , dtype = getNetworkDatatype_numpy() ))
    for i in range(len(interactionsToTest) +offset) :
        # produce an artifial person's SNP data with only the proposed interactions having values
        artificial_SNP_data =  np.zeros( (1,totalNumSNPs) , dtype = getNetworkDatatype_numpy())
        if i == len(interactionsToTest) : # if it is the last item, IE when all are interactions are active
            SNPlabels.append("All")
            for j in range(len(interactionsToTest) ) :
                SNP_set =  interactionsToTest[j]
                artificial_SNP_data[:,SNP_set] += strength ## add the SNPs in at each location
        else :
            SNP_set =  interactionsToTest[(i-1)] # all_possible_interactions[i]
            SNPlabels.append( np.array2string(SNP_set) )

        artificial_SNP_data[:,SNP_set] += strength

        targetLayerActivations = produceActivation(model,device,targetLayerIndex,artificial_SNP_data)
        if subtractNullActivations: targetLayerActivations -= nullActivations
        

        activationMap[i,:] = targetLayerActivations
           
    activationMap = np.abs(activationMap)
    # plot heatmap on grid: https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib   and    https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html
    os.makedirs(os.path.dirname(plotFolder_used), exist_ok=True)  
    fig_file = plotFolder_used +outFileStem

    fig, ax = plt.subplots(figsize=scale)
    if normalized is not None : ax.imshow(activationMap, cmap='gray', interpolation='nearest', vmin=normalized[0], vmax=normalized[1])
    else : ax.imshow(activationMap, cmap='gray', interpolation='nearest')
    plt.xlabel("Neurons")
    plt.ylabel("Interaction candidates")
    #plt.figure(figsize=(12,16))
    # ... and label them with the respective list entries
    ax.set_yticks(np.arange(len(SNPlabels)))
    ax.set_yticklabels(SNPlabels)
    ax.set_xticks(np.arange(len(neuronlabels)))
    ax.set_xticklabels(neuronlabels)

    titleTest = "Neuron activations for SNPs: " + outFileStem 
    if rSQ != -1 : titleTest += "( r^2: " + str( round(rSQ,2) ) + ")"
    ax.set_title(titleTest)
    
    fig.tight_layout()
    
    plt.show()
    fig.savefig(fig_file)  
    #fig_file = fig_file + "2"
    minValue = np.min(activationMap)
    maxValue = np.max(activationMap)
    normalized = [minValue, maxValue]
    return( normalized, activationMap )
#endregion
    
###############################################################################
# Hyperopt
#region Hyperopt
###############################################################################
minLayers=1
minConvLayers=1
minNeurons=20  #100
maxEpochs=100
minFilters=15
def parameter_space_pytorch():
    maxNeurons = 1000 #4000
    maxLayers = 10 # 20

        
    global minLayers ; global minNeurons ; global maxEpochs ; global minFilters ;
    # Set up a list for the parameter search space.
    space = {
    # the number of LSTM neurons
    'firstLayerSize': set_int_range('firstLayerSize', minNeurons,maxNeurons), #  +1 as otherwise Upper bond is exclusive, IE a range of 5-21 would only cover 5-20
    'epochs': set_int_range('epochs', 10,maxEpochs ), #  +1 as otherwise Upper bond is exclusive, IE a range of 5-21 would only cover 5-20  
    'hidCount': set_int_range('hidCount', minLayers,maxLayers), #  +1 as otherwise Upper bond is exclusive, IE a range of 5-21 would only cover 5-20
    'dropout': hp.uniform('dropout', 0.0, 0.9), 
    'learnRate': hp.uniform('learnRate', 0.000001, 0.01), #,  # loguniform
    'convLayers': set_int_range('convLayers',minConvLayers, 5),
    'convFilters': set_int_range('convFilters',minFilters, 75),
    #'optimizer': hp.choice('optimizer',[0,1, 2]), # which optimizer to use, 0 SGD, 1 ADAM and 2 AMSGrad   
    'hidAct': hp.choice('hidAct',[0,6])
      #'bnorm': hp.choice('bnorm',[0,1])
     # 'l2': hp.uniform('l2', 0.0, 20.0), # corresponds to h2 choice between 100% (delta = 0) to 5% (delta = 20)
    }
    return space


def set_int_range(name, myMin, myMax): # Set up a parameter range based on the given min and max values.
	# If the myMin and myMax values are equal, don't search over this parameter.
	if(myMin == myMax):
		return myMin

	# Swap the values so they are in the correct order if necessary.
	if(myMin > myMax):
		t=myMax
		myMax=myMin
		myMin=t
		
	# Randomly search over all integer values between the myMin and myMax
	return hp.choice(name, np.arange(myMin,myMax, dtype=int))
    

def plot_optimization_pytorch(trials, parameter, out_folder= ""): # , regression = False
    os.makedirs(os.path.dirname(out_folder), exist_ok=True) 
    # Create the base figure and axes.
    fig = plt.figure()
    ax = plt.subplot(111)

    # Pull out the trial data for the requested parameter.
    xs = [t['misc']['vals'][parameter] for t in trials.trials]
    ys = [-t['result']['loss'] for t in trials.trials]


    if parameter == 'convLayers':  # really hacky solution to start not at 0, but at the proper minimum
        #print("plotting layers, so we offset it")
        global minFilters
        for i in range(len(xs) ) :
            xs[i] = [s + minConvLayers for s in xs[i]] # element wise offset each value by the minimum value 


    if parameter == 'convFilters':  # really hacky solution to start not at 0, but at the proper minimum
        #print("plotting layers, so we offset it")
        global minFilters
        for i in range(len(xs) ) :
            xs[i] = [s + minFilters for s in xs[i]] # element wise offset each value by the minimum value 

    if parameter == 'hidCount':  # really hacky solution to start not at 0, but at the proper minimum
        #print("plotting layers, so we offset it")
        global minLayers
        for i in range(len(xs) ) :
            xs[i] = [s + minLayers for s in xs[i]] # element wise offset each value by the minimum value 
   
    if parameter == 'firstLayerSize':  # really hacky solution to start not at 0, but at the proper minimum
        #print("plotting layers, so we offset it")
        global minNeurons
        for i in range(len(xs) ) :
            xs[i] = [s + minNeurons for s in xs[i]] # element wise offset each value by the minimum value 

    # for nested/conditional params we need to remove values where parameter was disabled
    indicesToRemove = list()
    for i in range(len(xs) ) :
        if len(xs[i]) == 0 : 
            indicesToRemove.append(i) 
			#print("removing index ", i)

    # xs_2 = list( numpy.delete(xs, indicesToRemove) ) # this somehow 'flattens' nested lists ie [[0], [0]] becomes [0 , 0] NOT GOOD!
    xs = [i for j, i in enumerate(xs) if j not in indicesToRemove]
    ys = [i for j, i in enumerate(ys) if j not in indicesToRemove]
   
    for i in range( len(ys) ):
        if ys[i] < 0 : ys[i] = 0.0 # don't want to plot negative r^2s...

    #IDFileName=str(os.getcwd())+"/"+parameter+"_results.txt"
    IDFileName=out_folder+"/"+parameter+"_results.txt"
    with open(IDFileName, "w") as idFile: 
        idFile.write(parameter + "\t" + "response" + "\n")
        for i in range( len(xs) ):
            idFile.write( str(xs[i][0]) + "\t" + str(ys[i]) + "\n")

	# Draw the plot.
    ax.scatter(xs, ys, s=20, linewidth=0.01, alpha=0.5)
    ax.set_xlabel(parameter, fontsize=12)
    ylabel ='loss'
    #ylabel ='AUC'
    #if regression : ylabel = 'r^2'
    ax.set_ylabel(ylabel, fontsize=12)

    # Save the figure to file.
    fig_file=out_folder+"/"+parameter+"_optimisation.png"
    print("fig_file is: " + str(fig_file))
    fig.savefig(fig_file)


def mergeArgsAndParams(args,params) :
    argsMerged = copy.deepcopy(args) # createa deep copy of the original args object
    argsMerged = vars(argsMerged) # convert it to dictionary, so that we can easily copy the key/value pairs

    # go through the keys params and overwrite the corresponding entry in the args  (the argnames must match)
    for key, value in params.items():
        argsMerged[key] = value

    #argsMerged = { 'no_cuda': False, 'batch_size': 64, 'test_batch_size': 1000, 'epochs': 10, 'lr':0.01, 'momentum': 0.5, 'seed': 1, 'log_interval': 10 }
    argsMerged = SimpleNamespace(**argsMerged) # convert back to namespace, as that is what the build_model expects
    return(argsMerged)   
#params = {'epochs': 43, 'layers': 9, 'learnRate': 0.007538107949269826, 'neurons': 733, 'p_dropout': 0.7599124435352762}


numTrials_pytorch = 0
def trial_pytorch(params,device, args, y, M , train_minibatchIndices, PRS_betas , covars_all , test_minibatchIndices ):
    global multiGPU
    global numSNPs # cant use train_GWAS[0].shape[1], as in CNNs this will always be 1, as channels=1
    global supressOutput
    global numTrials_pytorch
    global numCovars
    numTrials_pytorch += 1

    # create model 
    argsMerged = mergeArgsAndParams(args,params) # produce a unified args object that has the 'on trial' parameters from hyperopt
    model = build_model(argsMerged, numSNPs, numCovars, y.shape[1],suppressPrint = True) # Don't pre-allocate memory; allocate as-needed ??
    model = modelToDevice(model,argsMerged,device, suppressPrint = True)

    # train model                        
    results = learn(model,device, argsMerged, y, M, train_minibatchIndices, PRS_betas=PRS_betas , covars_all =covars_all, test_minibatchIndices=test_minibatchIndices, saveModelLocation = None, epochMaxImproveThreshold = argsMerged.epochMaxImproveThreshold, learnRate =argsMerged.learnRate, momentum = argsMerged.momentum, suppressPrint = True, half = argsMerged.half, decayRate = argsMerged.LRdecay, accumulation_steps =argsMerged.accumulation_steps, l2= argsMerged.l2, optimizerChoice = argsMerged.optimizer, pheno_is_binary = pheno_is_binary)  #  (model,device, argsMerged, y, M, train_minibatchIndices, PRS_betas=PRS_betas , covars_all =covars_all, test_minibatchIndices=test_minibatchIndices, eval_train=True, eval_test=True, eval_freq = args.evalFreq, decayEnabled = False, suppressPrint = True)
    loss_trial = results['results']['valid_loss'][-1]
    lowestLoss=  results['results']['lowestLoss']
    lowestLoss_epoch=  results['results']['lowestLoss_epoch']

    if supressOutput == False: print('Trial ',str(numTrials_pytorch),' with parameters: ', str(params), " final loss: ", str(loss_trial), "but lowest loss was : " , str(lowestLoss), " at epoch: " , str(lowestLoss_epoch) )
    # params['epochs'] = int(lowestLoss_epoch) # overwrite this , this does NOT change the record in the trials.trails object

    if np.isnan(loss_trial) : 
        if supressOutput == False: print("loss is nan, set it to 0 ")
        loss_trial = 0
        
    attachments = {'lowestLoss':lowestLoss, 'lowestLoss_epoch': lowestLoss_epoch} 

    # Return the statistics for the best model in this trial.
    return {'loss': lowestLoss, 'status': STATUS_OK, 'attachments': attachments} 
  


def optimize_model_pytorch(device, args, y, M, train_minibatchIndices, PRS_betas = None, covars_all = None, test_minibatchIndices = None, out_folder ="", startupJobs = 40, maxevals = 200, noOut = False):
    global numTrials_pytorch

    numTrials_pytorch= 0

    trials = Trials()
    trial_wrapper = partial(trial_pytorch,device = device, args = args , y = y, M = M, train_minibatchIndices = train_minibatchIndices, PRS_betas = PRS_betas, covars_all = covars_all, test_minibatchIndices = test_minibatchIndices)

    best_pars = fmin(trial_wrapper, parameter_space_pytorch(), algo=partial(tpe.suggest, n_startup_jobs=(startupJobs) ), max_evals=maxevals, trials=trials)

    # Print the selected 'best' hyperparameters.
    if noOut == False: print('\nBest hyperparameter settings: ',space_eval(parameter_space_pytorch(), best_pars),'\n')

    # loops through the 1st entry in the dict that holds all the lookup keys
    #regression = True


    best_pars = space_eval(parameter_space_pytorch(), best_pars) # this turns the indices into the actual params into the valid aprameter space
    
    # override the epochs with the early start
    lowestLossIndex = np.argmin(trials.losses())
    trials.trial_attachments(trials.trials[lowestLossIndex])['lowestLoss_epoch']
    best_pars['earlyStopEpochs'] = trials.trial_attachments(trials.trials[lowestLossIndex])['lowestLoss_epoch']
    best_pars['earlyStopEpochs'] += 1 # as epochs are 0 based otherwise...
    best_pars['epochs'] = best_pars['earlyStopEpochs'] 
    if best_pars['epochs'] <= 0 : best_pars['epochs'] = 1 # we dont want a network without any training, as that will cause a problem for deep dreaming


    if noOut == False: print('\nBest hyperparameter settings (Final): ',best_pars,'\n')


    # do this last as this may crash the thing
    for p in trials.trials[0]['misc']['idxs']: plot_optimization_pytorch(trials, p, out_folder = out_folder) # , regression


    return(best_pars)


def writeKNeT_bestPars(outFile, best_pars ) :
    with open(outFile + ".best_pars", "w") as file: 
            file.write(str(best_pars) ) 


def loadKNeT_bestPars(outFile) :
    s=""
    with open(outFile + ".best_pars", "r") as file: 
        for i in file:
            s+=i
    best_pars= literal_eval(s)        
    return(best_pars) 
#endregion

###############################################################################
# Helper utils
#region Helper utils
###############################################################################  
# plots either a ROC or a Precision Recall
def plotROC_or_PrecisionRecall(out,xvals,yvals,auc_score, isROC = True, no_skill = 1.0, col = 'b', titleStem= "NN ") :
    if isROC :
        xaxis = 'False Positive Rate'
        yaxis = 'True Positive Rate'
        plotTitle = titleStem + "ROC"
    else :
        xaxis = 'Recall'
        yaxis = 'Precision'
        plotTitle = titleStem + "Precision-Recall"    
        
    fig = plt.figure()
    plt.title(plotTitle)
   
    plt.plot(xvals, yvals, col ,label= 'AUC= %0.3f'% (auc_score)) 
    noSkillLabel = 'No Skill= ' + '%0.3f'% (no_skill)
    if isROC : plt.plot([0, 1], [0, no_skill], linestyle='--', label=noSkillLabel)
    else : plt.plot([0, 1], [no_skill, no_skill], linestyle='--', label=noSkillLabel)

    plt.legend(loc='upper right')
    plt.xlim([-0.1,1.08])
    plt.ylim([-0.1,1.08])
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    plt.show()
    
    # Save the figure to file.
    fig_file=out+".pdf"
    print("written plot to:: " + str(fig_file))
    fig.savefig(fig_file)  

    fig_file=out+".png"
    fig.savefig(fig_file)  


def plotNNtraining(data,location, training_metric, valid_metric, acc_measure = "MSE", epoch_offset = 0) : # as the epoch's first training loss is much higher as that reflects before training accuracy, therefore we want to start at second epoch

    content = None
    cols = list()
    
    if platform.system().find('Windows') == -1 :
        print("Matlab uses Agg as we are on a *nix")
        matplotlib.use('Agg')
        
    import matplotlib.pyplot as plt
    
    if training_metric in data:
        #print("traing exists")
        cols.append('Traning')
        if(content is None) : content = data[training_metric][epoch_offset:len(data[training_metric])]
        else : content = np.column_stack( (content, data[training_metric][epoch_offset:len(data[training_metric])] ) )
        
    if valid_metric in data:
        #print("test exists")
        cols.append('Validation')
        if(content is None) : content = data[valid_metric][epoch_offset:len(data[valid_metric])]
        else : content = np.column_stack( (content, data[valid_metric][epoch_offset:len(data[valid_metric])]  ) ) 
               
    df = pd.DataFrame(content, index=data["epochs"][epoch_offset:len(data["epochs"])], columns=cols )

    plt.figure()
    
    ax = df.plot(title = "NN learning curve") 
    ax.set_xlabel("epochs")
    ax.set_ylabel(acc_measure)
    
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(location + '.eps', format='eps', dpi=1000)
    fig.savefig(location + '.png', dpi=300)



def printElapsedTime(start,end, text ="") : # https://stackoverflow.com/questions/27779677/how-to-format-elapsed-time-from-seconds-to-hours-minutes-seconds-and-milliseco
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print(text + "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds), flush=True)
        



def loadYMeanSds(location) : # attempts to load SD and mean for the y (text file assumed to have 2 lines, the mean and sd)
    y_mean = None
    y_sd = None
    
    my_file = Path(location)
    if my_file.is_file():
        with open(location, "r") as id:
            for i in id:
                itmp = i.rstrip().split()
                
                if y_mean is None : y_mean =  float(itmp[0])
                else : y_sd = float(itmp[0])

    return(y_mean,y_sd)
    
  
def writeYMeanSds(outFile, y_mean,y_sd ) :
    with open(outFile + ".ymeansd", "w") as file: 
            file.write(str(y_mean) + "\n" ) 
            file.write(str(y_sd) + "\n" ) 
 

def y_origscale(y,y_mean, y_sd): # transforms y back to raw scale
    y_rawscale = y *y_sd + y_mean 
    #np.mean(y_rawscale) # 2.0
    #np.var(y_rawscale) # 4.0
    return (y_rawscale)


def y_standardise(y,y_mean, y_sd):  # transforms a raw y into zscore
    y_standardised = (y - y_mean) / y_sd 
    #np.mean(y_standardised) # 2.0
    #np.var(y_standardised) # 4.0
    return (y_standardised)
#endregion

    # MORE CODE AND EXPLANATIONS FOR IN-PLACE NUMPY ARRAY REARRANGEMENT
    #region MORE CODE AND EXPLANATIONS FOR IN-PLACE NUMPY ARRAY REARRANGEMENT

    #PRS_SNPs= np.array(PRS_SNPs)
    #PRS_SNPs[ PRS_SNP_indices_in_bim ] = PRS_SNPs   #  PRS_SNPs_asArray = PRS_SNPs_asArray[PRS_SNP_indices_in_bim] this would be completely wrong!  PRS_SNPs_asArray[0] # [        2, 238368730], 
    ##PRS_betas[ PRS_SNP_indices_in_bim ] = PRS_betas
    #PRS_betas = PRS_betas[ PRS_SNP_indices_in_bim ] 




    #i=0
    #counter = 0
    #for i in range(len(bim_SNPs) ):
    #    equal =  PRS_SNPs[i][0] == bim_SNPs[i][0] and PRS_SNPs[i][1] == bim_SNPs[i][1]
    #    if(equal == False) :  
    #        counter+=1
    #    print(i, " PRS_SNPs: " , PRS_SNPs[i], " / bim_SNPs: " , bim_SNPs[i], "are equal: ", equal) 
            

    #print("not equal: " , counter, " / : " , len(bim_SNPs) )
    #   # not the same for 336 times wtf.....

    ## 336+485


    #counter = 0
    #for i in range(len(PRS_betas) ):
    #    equal =  PRS_betas2[i] == PRS_betas[i]
    #    if(equal == False) :  
    #        counter+=1
    #    print(i, " PRS_betas2: " , PRS_betas2[i], " / PRS_betas: " , PRS_betas[i], "are equal: ", equal) 
            

    #print("not equal: " , counter, " / : " , len(bim_SNPs) )










    #len(PRS_SNPs)
    #len(bim_SNPs)
    #for i in range(len(bim_SNPs) ):
    #    print("PRS_SNPs_asArray: " , PRS_SNPs[i], " / bim_SNPs: " , bim_SNPs[i]) # PRS_SNPs: ", PRS_SNPs[i], " / 


    #PRS_SNP_indices_in_bim[0] # 485
    #PRS_SNPs[0] # (10, 3820787) right..., so the first PRS SNP, is the 485th in the bim
    #PRS_SNPs
    #PRS_SNPs_matched = np.array(PRS_SNPs)[PRS_SNP_indices_in_bim]
    #PRS_SNPs_matched[0]
    #bim_SNPs[0]
    #PRS_SNPs[0]
    #len(PRS_SNPs)
    #bim_SNP_indices[0]
    #bim_SNPs[0] # (1, 2069172)
    #bim_SNPs[485] # (10, 3820787)
    #len(bim_SNPs) # 821
    #len(PRS_SNP_indices_in_bim)


    #counter = 0
    #for i in range(len(bim_SNPs) ):
    #    equal =  PRS_SNPs_asArray[i][0] == bim_SNPs[i][0] and PRS_SNPs_asArray[i][1] == bim_SNPs[i][1]
    #    if(equal == False) :  
    #        print(i, " PRS_SNPs: " , PRS_SNPs_asArray[i], " / bim_SNPs: " , bim_SNPs[i], "are not equal") 
    #        counter+=1

    #print("not equal: " , counter, " / : " , len(bim_SNPs) )



    #bim_SNPs[0] # the first SNP in the bim file is: (1, 2069172)
    #PRS_SNPs[0] # but the first SNP in the PRS file is: 10, 3820787)
    #PRS_SNP_indices_in_bim[0] # which is the 485th in the bim file: 485
    #bim_SNPs[485] # (10, 3820787)
    ## OK, now we want to move the 1st SNP in the PRS file to be the 485th, so that they are the same as the PLINK bim file
    #PRS_SNPs_asArray = np.array(PRS_SNPs)
    #PRS_SNPs_asArray[ PRS_SNP_indices_in_bim ] = PRS_SNPs   #  PRS_SNPs_asArray = PRS_SNPs_asArray[PRS_SNP_indices_in_bim] this would be completely wrong!  PRS_SNPs_asArray[0] # [        2, 238368730], 
    # PRS_SNPs[ PRS_SNP_indices_in_bim ] = PRS_SNPs # this would also be bad, as in-place rearranging numpy array doesnt work, I either need to use a custom function, or have an intermediary file. See answer #2: https://stackoverflow.com/questions/26239802/in-place-numpy-array-sorting-according-to-given-index
    
    # example with same number of elements
    #arr = np.array([10, 20, 30, 40, 50])
    #idx = [1, 0, 3, 4, 2]
    #arr[idx] # array([20, 10, 40, 50, 30])
    #arr2 = arr[idx] # this does work # array([20, 10, 40, 50, 30])
    #arr2
    #arr3 = np.array(arr)
    #arr3[idx] = arr # array([20, 10, 50, 30, 40])
    #arr3
    #arr4= np.array(arr)
    #arr4[idx] = arr[idx] # array([10, 20, 30, 40, 50]), does not work, as it just changes and restores original order
    #arr4
    #arr[idx] = arr # this does NOT work, as while we overwrite it we get duplicate entries: array([10, 10, 30, 30, 30])
    #arr
    ########################
    ## example with the idx being a subset
    #arr = np.array([10, 20, 30, 40, 50])
    #idx = [1, 0, 3, 4]
    #arr[idx] #array([20, 10, 40, 50])
    #arr2 = arr[idx] # this does work # array([20, 10, 40, 50])
    #arr2 # 
    #arr3 = np.array(arr)
    #arr3[idx] = arr # this does NOT work: valueError: shape mismatch: value array of shape (5,) could not be broadcast to indexing result of shape (4,)
    #arr3
    #arr4= np.array(arr)
    #arr4[idx] = arr[idx] # array([10, 20, 30, 40, 50]), this keeps the original, even when we have a subset, as we only change the subset indices, but keep the unchanged as is
    #arr4
    #arr5= np.zeros(len(idx),arr.dtype)
    #arr5[idx] = arr # this does NOT work: ValueError: shape mismatch: value array of shape (5,) could not be broadcast to indexing result of shape (4,)
    #arr5
    #arr[idx] = arr # this does NOT work: ValueError: shape mismatch: value array of shape (5,) could not be broadcast to indexing result of shape (4,)
    #arr
    ###########################################
    ########################

    ## example with the idx being a subset
    #arr = np.array(["SNP1", "SNP2", "SNP3", "SNP4", "SNP5"])
    ## find the indices of the SNPs in the PRS, now mapped to the subset of SNPs actually existing in both the PLINK and PRS files
    #dict_indices = {}
    #for i in range(len(arr) ): dict_indices[arr[i]] = i # save antoher lookup dict that stores the indices 


    #arr_sub = np.array(["SNP2", "SNP5", "SNP3"])
    #arr_sub_betas = np.array([0.3, 0.7, 0.1])
    #idx = matchIndices(arr_sub, dict_indices, True)
    #idx # array([1, 4, 2]) # this is the order according to the large array,
    #arr_sub[idx] # IndexError: index 4 is out of bounds for axis 0 with size 3

    #idx_v2 = matchIndices(arr_sub, dict_indices)
    #idx_v2 # array([0, 1, 2]) # this is the order according to the arr_sub, which isn't interesting
    ##  whereas, it should be 0,2,1, IE the indices referring to the arr_sub, but ordered according to the large array

    ## loop from the large array, check if element exists, and if yes, record its index
    #dict_indices_2 = {} # dictionary holds the lookups from the small array this time!
    #for i in range(len(arr_sub) ): dict_indices_2[arr_sub[i]] = i 
    #idx_v3 = matchIndices(arr, dict_indices_2, True)
    #idx_v3 # array([0, 2, 1]) # THIS IS CORRECT1!!
    #arr_sub[idx_v3] # array(['SNP2', 'SNP3', 'SNP5'], dtype='<U4') # reorder small array according to the large array,
    #arr_sub_betas[idx_v3] # rray([0.3, 0.1, 0.7]) # matching betas too!

    #matchIndices(arr, dict_indices_2) # array([1, 2, 4])
    # PRS_SNP_indices = matchIndices(PRS_SNPs, bim_SNPs_dict_indices)  # this would give me the match according to the PRS file... but I need to match it against the bim, so that I could reorder the PRS weights to match the PLINK order




    #i=0
    #counter = 0
    #for i in range(len(Indis_in_plink) ):
    #    equal =  Indis_in_plink [ i ]  == covarIDs_keeplist[i]
    #    if(equal == False) :  
    #        counter+=1
    #    print(i, " Indis_in_plink: " , Indis_in_plink [i ] , " / covarIDs_keeplist: " , covarIDs_keeplist[i], "are equal: ", equal) 
    #print("not equal: " , counter, " / : " , len(Indis_in_plink) )

    #i=0
    #counter = 0
    #for i in range(len(Indis_in_plink) ):
    #    equal =  Indis_in_plink [ i ]  == phenos_keeplist[i]
    #    if(equal == False) :  
    #        counter+=1
    #    print(i, " Indis_in_plink: " , Indis_in_plink [i ] , " / phenos_keeplist: " , phenos_keeplist[i], "are equal: ", equal) 
    #print("not equal: " , counter, " / : " , len(Indis_in_plink) )



    
    #np.random.seed(0)
    #dims = (1,5,5)
    #test_array = np.arange(np.prod(dims)).reshape(*dims)
    #test_array
    #test_array.shape
    #idx_dim1 = np.array([0,2,4])
    #idx_dim2 = np.array([1,3])
    #test_array[:,idx_dim1,idx_dim2]
    #test_array[:, idx_dim1[:,None], idx_dim2]

    

    #np.random.seed(0)
    #dims = (5,5)
    #test_array = np.arange(np.prod(dims)).reshape(*dims)
    #test_array
    #test_array.shape
    #idx_dim1 = np.array([0,2,4])
    #idx_dim2 = np.array([1,3])
    #test_array[idx_dim1,idx_dim2] # does not work
    #test_array[idx_dim1[:,None], idx_dim2] # works

    #test_array[np.ix_(idx_dim1,idx_dim2)] # works

    #test_array[list(idx_dim1), idx_dim2] # does not work



    #endregion

