# -*- coding: utf-8 -*-

#MIT License

#Copyright (c) 2019 Martin Kelemen

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
import numpy as np
import random
import copy
from itertools import combinations
import sys 
import time

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.nn.init as init
from torch.autograd import Variable

import operator
from collections import OrderedDict
from itertools import islice



###############################################################################
# Main vars
###############################################################################


OUT_MULTICLASS = "OUT_MULTICLASS"
OUT_REGRESSION = "OUT_REGRESSION"
OUT_MAE = "OUT_MAE"
EPSILON = 1e-8 # this is for FP32, for FP16, this sohuld be 1e-4
ALMOST_1 = 1.0 - EPSILON


OPTIMIZER_SGD = 0
OPTIMIZER_ADAM = 1
OPTIMIZER_AMSGRAD = 2

#device = None
NETWORK_DATATYPE = torch.float32# 'float32'
activations = list()
def getNetworkDatatype_numpy() :
    global NETWORK_DATATYPE
    if NETWORK_DATATYPE == torch.float32 : return np.float32
    else  : return np.float64

###############################################################################
# Learning
###############################################################################
#eval_train=True 
#eval_test=True 
#eval_freq = args.evalFreq, 
#decayEnabled = False
#decayEpochStart = 10
#suppressPrint = False
#gamma = 0.999
def learn(model, device, args      , y, M, train_minibatchIndices, PRS_betas=False     , covars_cont =False      , covars_factor =False        , test_minibatchIndices=False                , eval_train=False, eval_test=True, eval_freq = 100,  gamma = 0.999, decayEnabled = True, decayEpochStart = 10, suppressPrint = False):
    global NETWORK_DATATYPE
    torch.set_default_dtype(NETWORK_DATATYPE)
    model_training_orig = model.training
    setModelMode(model, True)
    
    #global device
  
    #if torch.cuda.is_available() and args.gpu ==1 : model.cuda() # mode model to cuda if available
    
    # I) setup optimiser & loss function
    if args.optimizer == OPTIMIZER_ADAM : optimizer = optim.Adam(model.parameters(), lr=args.learnRate, betas=(args.momentum, gamma), eps=EPSILON)   
    else : optimizer = optim.SGD(model.parameters(), lr=args.learnRate, momentum=args.momentum)
    # decide on loss function
    if y.shape[1] == 1 : criterion = nn.MSELoss() # loss function for regression is MSE, regression only 1 column
    else : criterion = nn.CrossEntropyLoss() 
    
    # cosmetics: if we are classifying then we are evaluating via 'accuracy' where as if we are regressing, then we care about error
    if  y.shape[1]  > 1 and isinstance(model[-1],nn.Softmax) : outPutType = OUT_MULTICLASS  # if its multicolumn AND we have softmax (IE  the cols sum to 1 for a prob), then its a multiclass classifcation problem
    elif  y.shape[1]  > 1 : outPutType = OUT_MAE # if it has more than 1 column, but the final layer isn't softmax, then we can only evaluate this by Mean Average Error
    else : outPutType = OUT_REGRESSION # if the thing we are trying to predict has only 1 column, then it is a regression problem            
    if outPutType == OUT_REGRESSION :  evaluation = "prediction (r^2)" 
    elif outPutType == OUT_MULTICLASS : evaluation ="accuracy (%)"  
    else : evaluation = "MAE"                 
    
    # setup results & early stop logic ( record highest validation accuracy and its epoch)
    results = {}
    results["epochs"] = list()
    results["train_accuracy"]  = list()
    results["test_accuracy"]  = list()
    results['highestAcc'] = -1.
    results['highestAcc_epoch'] = -1

    eta = args.learnRate
    decay = 0.0
    if decayEnabled : decay = 0.001 #eta / args.epochs
    if suppressPrint == False : print ("Pytorch start for #epochs: " + str(args.epochs) + " with LR decay enabled at: " + str(decay) , flush=True  )
    t = 0
    while t < args.epochs: #for t in range(0, args.epochs):
        out_str = " | it: "  + str(t) # "[{0:4d}] ".format(t)
        start_time = time.time()

        # 1) Complete an entire training cycle: Forward then Backward propagation, then update weights, do this for ALL minibatches in sequence
        currentBatchNum = 0
        batchNum = 0
        totalBatches = len(train_minibatchIndices)
        for batchNum in range(len(train_minibatchIndices)):
            b_labels, b_data = getMiniBatchData(y, M, train_minibatchIndices[batchNum], PRS_betas,  covars_cont, covars_factor)
            #b_data = train_X[batchNum] ; b_labels = train_Y[batchNum]
            #print ("batch: " + str(currentBatchNum) , flush=True  )
           
            # convert data to torch & move it to CUDA if available
            b_data = dataToDevice(b_data, device)
            b_labels = torch.from_numpy(b_labels).to(device)

            # perform full learning cycle, FP, BP and update weights
            #learnCycle(model, criterion, optimizer,args, b_data, b_labels)
            optimizer.zero_grad()   # zero the gradient buffers
            resetCachedActivations() # reset any saved activations before FP
            
            # Forward Propagate
            yhat = model(b_data) 
            loss = criterion(yhat, b_labels) # loss function
            
            # apply regularization to loss
            addRegularizerCosts(model, loss, args)
        
            # Backprop
            loss.backward()  
            
            # update Weights
            optimizer.step() 
            
            # update prograss bar
            barPos =  currentBatchNum / totalBatches # the bar position in %
            barPos = round(20 * barPos) # scale it to 1-20
            if suppressPrint == False : 
                sys.stdout.write('\r')
                sys.stdout.write("[%-20s] %d%%" % ('='*barPos, 5*barPos))
                sys.stdout.flush()  
            currentBatchNum += 1.0
        
 
        if t % eval_freq == 0 or t == args.epochs -1: # only evaluation fit every 100th or so iteration, as that is expensive
            results["epochs"].append(t) # add the epoch's number
            
            if eval_train or t == args.epochs -1:
                accuracy = evalAccuracy(model, device, args, outPutType, y, M, train_minibatchIndices, PRS_betas,  covars_cont, covars_factor) # evalAccuracy(model, device, args, outPutType, train_X, train_Y)  
                out_str =  out_str + " / Training " +evaluation +": "+ str( accuracy )                  
                results["train_accuracy"].append(accuracy)  

            if test_minibatchIndices is not None: # if a test set was supplied at all
                if eval_test or t == args.epochs -1:
                    accuracy = evalAccuracy(model, device, args, outPutType, y, M, test_minibatchIndices, PRS_betas,  covars_cont, covars_factor) # evalAccuracy(model, device, args, outPutType, test_X, test_Y, validEval = True) 
                    
                    if accuracy > results['highestAcc'] :
                        results['highestAcc'] = accuracy
                        results['highestAcc_epoch'] = t
                        
                    results["test_accuracy"].append(accuracy)
                    out_str =  out_str + " / Test " +evaluation +": "+  str( accuracy )
      
        #gc.collect()
        elapsed_time = time.time() - start_time 
        if suppressPrint == False : print(out_str + " / " + str( round(elapsed_time) ) + " secs (LR: " + str(eta) + ")" , flush=True)
        
        # update learning rate
        if t > decayEpochStart : eta = eta * 1/(1 + decay * t) # as t starts as 0, this will only kick in for the 3rd iteration

        t += 1
    setModelMode(model, model_training_orig)
    return ( { "results" : results})  
 


###############################################################################
#Learning helpers
#region Learning helpers
###############################################################################


def addL2Regularizer(loss,model, regParam, targetLayers) : # add L2 to the first layer: https://discuss.pytorch.org/t/how-does-one-implement-weight-regularization-l1-or-l2-manually-without-optimum/7951/12
    for i in range( len(targetLayers) ) :
        layerParams =   model[targetLayers[i]].named_parameters() 
        for W in layerParams:  # dont regularize bias params  
            if 'bias' not in W[0]: loss += 0.5 * regParam * torch.pow(W[1], 2).sum()  
      
     
def addL1Regularizer(loss,model, regParam, targetLayers) :
    for i in range( len(targetLayers) ) :
        layerParams =   model[targetLayers[i]].named_parameters() 
        for W in layerParams:  # dont regularize bias params
            if 'bias' not in W[0]: loss +=  0.5 * regParam * torch.abs(W[1]).sum()   


                

def get_activation():
    global activations
    def hook(model, input, output):
        if model.training : # ONLY add this if the model is in training mode, bad idea, as we only 
            activations.append(output)  # output.detach()
    return hook



        
def resetCachedActivations() : # must be called BEFORE FP
    global activations
    activations = list() # reset this


        
    
def addRegularizerCosts(model, loss, args, subsetModelLayers = None) : # subsetModelLayers is a LIST of layers
    # findPrevNextWeightedLayer gets into an infinite loop for the tree NN, so we skip this 
    if args.hidl2 > 0.0 or args.l1 > 0.0 :
        penultimateIndex = findPrevNextWeightedLayer(model, startIndex = -1, startIndexOK = False)    # we want the first Linear layer with the Weights, starting from the output
    
        if args.hidl2 > 0.0 :  # if there are no subset layers, or there ARE but the actual layer exist in that list
            # for tree NNs the L2 gets applied to the first FC layer, which is NOT index 0
            firstFCLayerIndex = getFirstFCLayerIndex(model)
            if subsetModelLayers is None or subsetModelLayers is not None and model[firstFCLayerIndex] in subsetModelLayers : addL2Regularizer(loss, model, args.hidl2, [firstFCLayerIndex]) # L2 norm gets added to the first layer only      
        if args.l1 > 0.0 : 
            if subsetModelLayers is None or subsetModelLayers is not None and model[penultimateIndex] in subsetModelLayers : addL1Regularizer(loss, model, args.l1, [penultimateIndex]) # this gets added to the penultimate (inference) layer

#endregion
   
###############################################################################
#Inference functions:
#region Inference functions:
###############################################################################
def findInteractionsExhaustive(model, device, totalNumSNPs,  strength = 1, orderInteraction= 3) :
    global NETWORK_DATATYPE
    allInteractions = {}
    possible_SNPs = np.asarray( range(totalNumSNPs) ) 
    all_possible_interactions = combinations(possible_SNPs, orderInteraction) # np.array( list(  combinations(possible_SNPs, orderInteraction) ) )
    allStrengths = list() # np.zeros(len(all_possible_interactions))
    model_training_orig = model.training
    setModelMode(model, False)
   
    for i in all_possible_interactions : #   print(i)
        SNP_set =   np.array(i) # all_possible_interactions[i]
        # produce an artifial person's SNP data with only the proposed interactions having values
        artificial_SNP_data = np.zeros(totalNumSNPs, dtype = getNetworkDatatype_numpy())
        artificial_SNP_data[SNP_set] = strength
        artificial_SNP_data = artificial_SNP_data[np.newaxis,] # make sure its 2D
        artificial_SNP_data = dataToDevice(artificial_SNP_data, device)

        # produce phenotype prediction for this guy
        interactionStrength = model(artificial_SNP_data).detach().cpu().numpy()[0,0]
  
        allStrengths.append(interactionStrength)
        interactionStrength = np.abs(interactionStrength)
        SNP_set.sort() # enforce same order, so we can aggregate strengths for same interactions
        SNP_set = tuple(SNP_set) # need to convert to tuple otherwise cant use it as key for a dictionary
        allInteractions[SNP_set] = interactionStrength

    allInteractions_sorted = sorted(allInteractions.items(), key=operator.itemgetter(1))
    allInteractions_sorted.reverse() # start with the highest strength
    setModelMode(model, model_training_orig)
    allStrengths = np.array( allStrengths )
    #plt.hist(allStrengths, bins = 'auto', density = True)  # , 200)
    return(allInteractions_sorted, allStrengths)


 # Deep Dreaming via gradient ascent: https://github.com/shksa/My-cs231n-assignment/blob/master/ImageGradients.ipynb   / https://github.com/XavierLinNow/deepdream_pytorch/blob/master/deepdream.py
def dream(model, device,args, target_class, targetClass_strength = 1, startGenome = None, iterations = 2,lr = 1.5 , startLayerIndex = 0, finishLayerIndex = -1, suppressPrint = True):
    model = getModel(model)
    if type(model) is not nn.Sequential : raise ValueError("DREAMING DOES NOT WORK FOR TREE NN, as b_data is a list here, also have to be mindful of how the module list has the regions as layers, IE with a tree of 2 regions, both activated will be stored as [Linear,Relu,Linear,Relu]  ...")
    #global device
    global NETWORK_DATATYPE
    
    # determine the shape of the synthetic genome to be generated, if startLayerIndex ==0, then, this will be the shape of the input, otherwise, it needs to be shaped as the input expected by that start layer
    if type(model[startLayerIndex]) == nn.Conv1d :  INPUT_DIM = [model[startLayerIndex].in_channels,model[startLayerIndex].out_channels ]  
    else :  INPUT_DIM = [model[startLayerIndex].in_features] # this is the number of SNPs
    
    # if we want to start from 0s or a supplied genome
    if startGenome is None : X_fooling = np.zeros( (1,*INPUT_DIM)  , dtype=getNetworkDatatype_numpy())  # start from blackness
    else: X_fooling = startGenome.copy() # start from a supplied genome


    # generate true labels: depending on the finishLayer's type:
    finishLayerIndex_used = findPrevNextWeightedLayer(model, finishLayerIndex) # need another variable to store this, so we can keep the orig as a flag if 'we are in the last layer mode'
    if type(model[finishLayerIndex_used]) == nn.Conv1d :
        labels = np.zeros( ( 1,*getConv1DOutput(model[finishLayerIndex_used]) ) , dtype=getNetworkDatatype_numpy())
        labels[:,target_class,:] = targetClass_strength # index [1] is the index for neurons 
    else: # if the final layer is a linear (this cannot be a softmax anyway as we are interested in the score not the actual classification probability)
        labels = np.zeros( (1, model[finishLayerIndex_used].out_features) , dtype=getNetworkDatatype_numpy())
        labels[:, target_class] = targetClass_strength

        
    mask_indices = np.where(labels == 0) # get indices where we would zero out activations (IE all but the targetClass),


    # convert data to torch & move it to CUDA if available
    b_data = dataToDevice(X_fooling, device)
    b_data.requires_grad = True # this is so that the Input also has Derivatives, otherwise BP would stop at the 1st layer
    b_labels = torch.from_numpy(labels).to(device)

    # setup hook for BP  
    #errorHook = model[startLayerIndex].register_backward_hook(...)  # add/remove hook: https://github.com/pytorch/pytorch/issues/5037
    
    # FP / BP only from/to the start/finish ( we want to FP/BP to the activation though as the crossentropy loss epxect probs, not scores)
    origLayers = list(model)
    lastLayerSlice = finishLayerIndex_used +1 # finds the nearest previous Linear layer, this could be last layer if it was a regression, +1 as this is a slice, that is EXCLUSIVE, IE 0:len
    subsetModelLayers = origLayers[startLayerIndex:lastLayerSlice]
    modelInference = nn.Sequential(*subsetModelLayers) # create a subset model of only
    modelInference.eval() # otherwise dropout and other layers would act up

    criterion = dream_loss() 
    #subsetModelLayers
    for i in range(iterations) :     
        if suppressPrint == False : print("dreaming iteration: "  + str(i))
 
       # learnCycle(modelInference, criterion, None,args, b_data, b_labels, updateWeights = False) # dont change weights

        # Forward Propagate
        resetCachedActivations() # reset any saved activations before FP, this flushes any activations
        yhat = modelInference(b_data) 
        if finishLayerIndex != -1 : yhat[mask_indices] = 0.0 #  this is only applied when we are interested in a particular neuron in NOT the final layer 
  
        loss = criterion(yhat, b_labels) # loss function
        
        # apply regularization to loss
        addRegularizerCosts(model, loss, args, subsetModelLayers = subsetModelLayers) # we must use the original 'model' and not the 'modelInference' here, as the latter may have the indices referring to different layers than what we used theregularizer for
    
        # Backprop
        loss.backward()          

        # ge the Error at the input
        dX= b_data.grad.data.detach().cpu().numpy().copy() # the error is saved directly onto this torch variable,  and the shape of the derivative is the SAME as the input
        b_data = b_data.detach().cpu().numpy().copy() # convert this to numpy
        
        # basic update, we keep creating the synthetic input which results in the highest activation
        b_data += lr * dX

        # convert back to torch for the next iteration
        b_data = dataToDevice(b_data, device)
        b_data.requires_grad = True # this is so that the Input also has Derivatives, otherwise BP would stop at the 1st layer
        
    # remove hook
    #errorHook.remove()
    
    modelInference.train() # reset model into traning mode
    return(b_data.detach().cpu().numpy().ravel())

    
#  evaluate Knet via associating each SNP via the Garson weights produces 'expected' importance scores from the Weight matrices of a network: this is ~ to the Importance scores from the 'deep dreaming' (which are the observed) but seem to be less accurate
def dream_Garson(model , startLayerIndex = 0, finishLayerIndex = -1 , activate = False) :
    model = getModel(model)
    if type(model) is not nn.Sequential : raise ValueError("DREAMING DOES NOT WORK FOR TREE NN, as b_data is a list here, also have to be mindful of how the module list has the regions as layers, IE with a tree of 2 regions, both activated will be stored as [Linear,Relu,Linear,Relu]  ...")
    
    model_training_orig = model.training
    setModelMode(model, False)

    importanceScores = None
    finishLayerIndex_used = findPrevNextWeightedLayer(model, finishLayerIndex)
    for layerIndex in range(startLayerIndex, (finishLayerIndex_used +1) ) :

        if type(model[layerIndex]) == nn.Conv1d or type(model[layerIndex]) == nn.Linear  or type(model[layerIndex]) == nn.Conv2d :
            modelParams = list(model[layerIndex].parameters() )
            weights = modelParams[0].t().detach()  # https://discuss.pytorch.org/t/clone-and-detach-in-v0-4-0/16861/2
           
            if activate :  # if we meant to activate weights, and there is a next layer which is an activation type of layer
                nextActivationIndex =  findPrevNextActivationLayer(model,layerIndex, prev = False, startIndexOK = False)
                if nextActivationIndex != -1: weights = model[nextActivationIndex](weights)
                    
            weights = weights.numpy()
            weights = np.abs(weights) # only take abs AFTER the nonlinearity
            if importanceScores is None :   importanceScores = weights
            else : importanceScores = importanceScores @ weights

    setModelMode(model, model_training_orig)    
    return(importanceScores.ravel())


    
# use the basic yhat-y loss: and not the l2 norm https://github.com/XavierLinNow/deepdream_pytorch/blob/master/deepdream.py
class dream_loss(torch.nn.Module): # custom loss: https://spandan-madan.github.io/A-Collection-of-important-tasks-in-pytorch/     ,,    https://discuss.pytorch.org/t/build-your-own-loss-function-in-pytorch/235/42
    def __init__(self):
        super(dream_loss,self).__init__()
        
    def forward(self, yhat, y):
        diff = torch.sum( (yhat - y) )
        return diff

#endregion
###############################################################################
# Helper utils
#region # Helper utils
###############################################################################

# Conv1D definition [in_channels, out_channels, kernel_size]  # in_channels = # SNPs, out_channels = number of neurons/filters
# Conv1D expected input shape [batch-size, in_channels, out_channels] # 
def getConv1DOutput(myCov1D) : # get the shape
    Cout = myCov1D.out_channels
    Lout = int ( ( myCov1D.in_channels + 2*myCov1D.padding[0] - myCov1D.dilation[0]*(myCov1D.kernel_size[0] -1) -1 ) / myCov1D.stride[0] +1 ) # https://pytorch.org/docs/stable/nn.html#torch.nn.Conv1d
    return( [Cout,Lout] )

# WARNING: when going FORWARD, these 2 functions (findPrevNextWeightedLayer, findPrevNextActivationLayer) may not work as expected for a tree type of NN, as the next layer in the module list could be a 'subregion' of a level (IE we have 2 regions, both with activations, these are stored as Linear,Relu,Linear,Relu, so we will find the Relu after the first region, and not the next 'level's activation)
    # but they should work OK if we want to find the 'last' layer as the output of it is the same
def findPrevNextWeightedLayer(model,startIndex = -1, prev = True, startIndexOK = True): # finds the next/prev 'proper' layer with weights
    model = getModel(model)
    layersAsList = list(model)
    if startIndex == -1 : startIndex = len(layersAsList) -1 # need to find the actual layer index if we were passed in the python shortcut of -1 
    
    if prev : step = -1 # depending if we want the next or previous layer, the step will be +1 or -1 to forward or backwards
    else : step = 1

    currentIndex = startIndex 
    if startIndexOK == False : currentIndex += step # if we can use the start index, then start from there
    while True:
        if currentIndex >= len(layersAsList) or currentIndex < 0: 
            raise Exception("run out of layers!")
            break
        if type(model[currentIndex]) == nn.Conv1d or type(model[currentIndex]) == nn.Linear  or type(model[currentIndex]) == nn.Conv2d : 
           # print("found layer at: ", currentIndex)
            break
        currentIndex += step
        print("currentIndex, is:" , currentIndex, " (out of num layers:", len(layersAsList) ,")for type(model[currentIndex]):" , type(model[currentIndex]))

    return(currentIndex)


def isLayerActivation(layer) :
    if type(layer) == nn.Sigmoid or type(layer) == nn.ReLU or type(layer) == nn.LeakyReLU or type(layer) == nn.Softplus or type(layer) == nn.SELU : return(True)
    else : return(False)

def findPrevNextActivationLayer(model,layerIndex, prev = True, startIndexOK = True) : # finds the next activation type layer (RELU/leaky relu / SELu etc)
    if prev : step = -1 # depending if we want the next or previous layer, the step will be +1 or -1 to forward or backwards
    else : step = 1
    model = getModel(model)
    layersAsList = list(model)
    currentIndex = layerIndex 
    if startIndexOK == False : currentIndex += step # if we can use the start index, then start from there
    while True:
        if currentIndex >= len(layersAsList) or currentIndex < 0:  
            print("No Activation found!")
            currentIndex = -1
            break
        if (isLayerActivation(model[currentIndex])) : break
        #if type(model[currentIndex]) == nn.Sigmoid or type(model[currentIndex]) == nn.ReLU or type(model[currentIndex]) == nn.LeakyReLU or type(model[currentIndex]) == nn.Softplus or type(model[currentIndex]) == nn.SELU : 
           # print("found layer at: ", currentIndex)
        #    break 
        currentIndex += step
        
    return(currentIndex)


def setModelMode(model, training = True): # just setting model.training = False, is NOT enough to set all layers to be in 'eval' mode, it will only set the wrapper
    if training : model.train()
    else : model.eval()


#batchyIndices = train_minibatchIndices[0]
def getMiniBatchData(y, M, batchIndices, PRS_betas,  covars_cont, covars_factor) :
    b_data = M[batchIndices]
    
    # multiply out the SNPs with the PRS weights (if any) at float32 precision
    if PRS_betas is not None : 
        b_data = b_data * PRS_betas #   b_data @ PRS_betas: this would NOT be what we want, as this would create a vector of size (minibatch,1), ie the per indi PRS. But here we just want each SNP [0,1,2] to be multiplied by their respective GWAS weight. sum(b_data[0]) would be then recover the PRS  b_data2[0] ==  sum(b_data3[0])
    else : b_data = b_data.astype('float32') # still need to cast this, in case we dont have either PRS betas or covars

    # depending on if we have covars, we concat that into the data
    if covars_cont is not None: 
        b_data_all = np.zeros( (b_data.shape[0],b_data.shape[1] + covars_cont.shape[1] ), dtype = np.float32 ) # create blank 2D array that could accommodate both the SNP and covars data
        b_data_all[:,0:b_data.shape[1]] = b_data # paste in the SNP data in the first part
        b_data_all[:,b_data.shape[1]:b_data.shape[1] + covars_cont.shape[1]] = covars_cont[batchIndices] # paste in the covars data in the second half
        b_data = b_data_all
        del b_data_all; gc.collect();

    # get the labels
    b_labels = y[batchIndices]

    # TODO: deal with covars_factor[batchyIndices]

    return b_labels, b_data


# common function that returns metric on how well a model does


def evalAccuracy(model, device, args, outPutType, y, M, minibatchIndices, PRS_betas,  covars_cont, covars_factor) :  # evalAccuracy(model, device, args, outPutType, test_X, test_Y, validEval = False) :
    model_training_orig = model.training
    setModelMode(model, False)
    N_test = len(minibatchIndices)*len(y[0]) # len(test_Y)*len(test_Y[0])
    totals = 0

    for batchNum in range(len(minibatchIndices)):
    #for b_data, b_labels in zip(test_X, test_Y):
        b_labels, b_data = getMiniBatchData(y, M, minibatchIndices[batchNum], PRS_betas,  covars_cont, covars_factor)
        b_data = dataToDevice(b_data, device) 
        b_labels = torch.from_numpy(b_labels).to(device) 
        b_labels = b_labels.view(b_labels.shape[0],-1 )  # make it the same shape as output
       
        yhat = model(b_data) # need to compute the Yhat again, as this is the yhat AFTER updating the weights, not before as in 'learn()' function
        b_data = None
        
        # depending on if we are in a classification or regression problem, we evaluate performance differently
        if outPutType == OUT_REGRESSION :  
            currentRate = torch_pearsonr( yhat.view(-1)  , b_labels.view(-1))**2   
            N_test = len(minibatchIndices) # as we are testing correlation, the N should refer to the number of batches, and NOT the total number of observations  
        elif outPutType == OUT_MULTICLASS : 
            currentRate = calc_Accuracy(yhat,b_labels )    # calculate accuracy, this is ROUNDED 
        else : # mean absolute error
            currentRate = -torch.mean(torch.abs(yhat - b_labels)) # negative as the rest of the metrics are accuracy, IE the greater the error, the lower the accuracy
            N_test = len(minibatchIndices) # as we are testing correlation, the N should refer to the number of batches, and NOT the total number of observations
 
        currentRate = float(currentRate.detach().cpu().numpy() )  # need to move it back to CPU
        totals = totals +currentRate # sum in all minibatches

    accuracy = round( float(totals)/N_test,5)
    setModelMode(model, model_training_orig)
    return(accuracy)

def calc_Accuracy(yhat, y): # calculates how many times we got
    y_preds= torch.argmax(yhat, axis=1)
    y_trues= torch.argmax(y, axis=1)
    num_matches = torch.sum(y_preds == y_trues)
    return(num_matches) 
    
    
def torch_pearsonr(x, y):  # https://github.com/pytorch/pytorch/issues/1254
    mean_x = torch.mean(x)
    mean_y = torch.mean(y)
    xm = x.sub(mean_x)
    ym = y.sub(mean_y)
    r_num = xm.dot(ym)
    r_den = torch.norm(xm, 2) * torch.norm(ym, 2)
    r_val = r_num / r_den
    return r_val

   
# using kaiming instead of xavier as it is better for RELUs: https://stackoverflow.com/questions/48641192/xavier-and-he-normal-initialization-difference
def weight_init(m): # https://gist.github.com/jeasinema/ed9236ce743c8efaf30fa2ff732749f5
    if isinstance(m, nn.Conv1d):
        init.normal_(m.weight.data)
        if m.bias is not None:
            init.normal_(m.bias.data)
    elif isinstance(m, nn.Conv2d):
        init.kaiming_normal_(m.weight.data)
        if m.bias is not None:
            init.normal_(m.bias.data)
    elif isinstance(m, nn.Conv3d):
        init.kaiming_normal_(m.weight.data)
        if m.bias is not None:
            init.normal_(m.bias.data)
    elif isinstance(m, nn.ConvTranspose1d):
        init.normal_(m.weight.data)
        if m.bias is not None:
            init.normal_(m.bias.data)
    elif isinstance(m, nn.ConvTranspose2d):
        init.kaiming_normal_(m.weight.data)
        if m.bias is not None:
            init.normal_(m.bias.data)
    elif isinstance(m, nn.ConvTranspose3d):
        init.kaiming_normal_(m.weight.data)
        if m.bias is not None:
            init.normal_(m.bias.data)
    elif isinstance(m, nn.BatchNorm1d):
        init.normal_(m.weight.data, mean=1, std=0.02)
        init.constant_(m.bias.data, 0)
    elif isinstance(m, nn.BatchNorm2d):
        init.normal_(m.weight.data, mean=1, std=0.02)
        init.constant_(m.bias.data, 0)
    elif isinstance(m, nn.BatchNorm3d):
        init.normal_(m.weight.data, mean=1, std=0.02)
        init.constant_(m.bias.data, 0)
    elif isinstance(m, nn.Linear):
        init.kaiming_normal_(m.weight.data)
        init.normal_(m.bias.data)
    elif isinstance(m, nn.LSTM):
        for param in m.parameters():
            if len(param.shape) >= 2:
                init.orthogonal_(param.data)
            else:
                init.normal_(param.data)
    elif isinstance(m, nn.LSTMCell):
        for param in m.parameters():
            if len(param.shape) >= 2:
                init.orthogonal_(param.data)
            else:
                init.normal_(param.data)
    elif isinstance(m, nn.GRU):
        for param in m.parameters():
            if len(param.shape) >= 2:
                init.orthogonal_(param.data)
            else:
                init.normal_(param.data)
    elif isinstance(m, nn.GRUCell):
        for param in m.parameters():
            if len(param.shape) >= 2:
                init.orthogonal_(param.data)
            else:
                init.normal_(param.data) 

        
def getModel(model) : # in case we need to access the underlying model of a DatParallelised model, dont know if this will return a model that has all the Weights correct or not...
    if type(model ) == nn.DataParallel : return model.module
    else: return model


def getFirstFCLayerIndex(model) : # in case we need to access the underlying model of a DatParallelised model, dont know if this will return a model that has all the Weights correct or not...
    model = getModel(model)
    if type(model ) == nn.Sequential : return(0) # for sequential models the first usable layer where we can apply h2, is always the 0th index
    else : return (model.firstFCLayerIndex)  # if it not sequential then it must be a 'Tree model', then we have saved this into a variable
        
#import sys

def dataToDevice(b_data, device) : # cast data to the device, if this is a 'tree' type NN, then the data will be a list of arrays, so we cannot directly cast it to cude all at once, but have to loop through and do it one by one
    b_data =  torch.from_numpy(b_data).to(device) 
    
    return(b_data)

#endregion