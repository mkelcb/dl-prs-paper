# -*- coding: utf-8 -*-

#MIT License

#Copyright (c) 2023 Martin Kelemen

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


#OUT_MULTICLASS = "OUT_MULTICLASS"
#OUT_REGRESSION = "OUT_REGRESSION"
#OUT_MAE = "OUT_MAE"
#EPSILON = 1e-8 # this is for FP32, for FP16, this sohuld be 1e-4
#ALMOST_1 = 1.0 - EPSILON


OPTIMIZER_SGD = 0
OPTIMIZER_ADAM = 1
OPTIMIZER_AMSGRAD = 2

#device = None
#NETWORK_DATATYPE = torch.float32# 'float32'
activations = list()


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


           # the sequence classifier:
         #(model, device, trainLoader, validLoader                                                                                                                                           , saveModelLocation = None, epochMaxImproveThreshold = 30, learnRate = 0.1, momentum = 0.9,  gamma = 0.999, suppressPrint = False, half = False, decayRate = 0.96, accumulation_steps =6, debugOut = None, l2= 0., sgd = False, pheno_is_binary = True, t = 0, training_hasntmprovedNumIts = 0, results = None, optimizerState = None):
         #(model, device, args      , y, M, train_minibatchIndices, PRS_betas=False     , covars_cont =False      , covars_factor =False        , test_minibatchIndices=False                , eval_train=False        , eval_test=True               , eval_freq = 100,  gamma = 0.999, decayEnabled = True, decayEpochStart = 10, suppressPrint = False):
         # the old knet       

#model
#device
#args
#y
#train_minibatchIndices
#M
#PRS_betas=PRS_betas 
#covars_all =covars_all
#test_minibatchIndices=test_minibatchIndices
#saveModelLocation = saveModelLocation

#epochMaxImproveThreshold = args.epochMaxImproveThreshold
#learnRate =args.learnRate
#momentum = args.momentum
#half = args.half
#decayRate = args.LRdecay
#accumulation_steps =args.accumulation_steps
#l2= args.l2
#optimizerChoice = args.optimizer
#pheno_is_binary = pheno_is_binary
#t = 0
#training_hasntmprovedNumIts = 0
#results = None
#optimizerState = None


def learn(model, device, args      , y, train_minibatchIndices, M=None, PRS_betas=None     , covars_all =None       , test_minibatchIndices=None                , saveModelLocation = None, epochMaxImproveThreshold = 30, learnRate = 0.1, momentum = 0.9,  gamma = 0.999, suppressPrint = False, half = False, decayRate = 0.96, accumulation_steps =6, debugOut = None, l2= 0., optimizerChoice = 0, pheno_is_binary = True, t = 0, training_hasntmprovedNumIts = 0, results = None, optimizerState = None) :
    model_training_orig = model.training
    setModelMode(model, True)
    numtrainingBatches =  len(train_minibatchIndices)

    # I) setup optimiser & loss function
    if half : EPSILON = 1e-4 # for FP16, this sohuld be 1e-4
    else : EPSILON = 1e-7  # 1e-8 default for pytorch,  1e-7 , default for keras

    if optimizerChoice ==  OPTIMIZER_SGD : optimizer = optim.SGD(model.parameters(), lr=learnRate, momentum=momentum, weight_decay= l2)
    else : optimizer = optim.Adam(model.parameters(), lr=learnRate, betas=(momentum, gamma), eps=EPSILON, weight_decay= l2)   
    
    if optimizerState is not None : # if resuming model, need to restore the optimizer too
        print("resuming optimizer")
        optimizer.load_state_dict(optimizerState)

    # decide on loss function
    if pheno_is_binary : criterion = nn.BCELoss() 
    else :  criterion = nn.MSELoss() # nn.MSELoss()

    # setup results & early stop logic (record highest validation accuracy and its epoch)
    if results is None :
        results = {}
        results["epochs"] = list()
        results["train_loss"]  = list()
        results["valid_loss"]  = list()
        results['lowestLoss'] = None
        results['lowestLoss_epoch'] = -1

    eta = learnRate
    if suppressPrint == False : print ("DL GXE PRS (EPSILON = ", EPSILON,"): Training started iterations until no improvement epochs: " + str(epochMaxImproveThreshold), "/ learn Rate:", str(eta), " / exponential decay:", decayRate, " / l2: ", l2, " / pheno_is_binary: ", pheno_is_binary, flush=True  )

    if decayRate != -1 : lr_scheduler = optim.lr_scheduler.ExponentialLR(optimizer=optimizer, gamma=decayRate)

    # if we are resuming a model and starting from a non 0 epoch, we need to 'catch up' the learning rate scheduler
    if decayRate != -1 and t > 0 :
        for i in range(t -1) :
            lr_scheduler.step()
            print(i,": resuming learning rate", lr_scheduler.get_last_lr()[0])

    while training_hasntmprovedNumIts < epochMaxImproveThreshold : # keep training until the model stops improving
        out_str = " | epoch: "  + str(t) # "[{0:4d}] ".format(t)
        results["epochs"].append(t) # add the epoch's number
        
        # 1) Complete an entire training cycle: Forward then Backward propagation, then update weights, do this for ALL minibatches in sequence
        start_time = time.time()
        currentBatchNum = 0
        epoch_trainingLoss = 0.   
        torch.set_grad_enabled(True) ; model.train() # enable grads for training
        if t == 0 : 
            torch.set_grad_enabled(False)
            model.eval() # for the first epoch we want to know what the baseline prediction looks like, must use .eval() otherwise BN-like layers would still learn/update their params which could cause numerical issues
        
       # yhats_saved = [] ;  b_labels_saved = []
       #batchNum=0
        optimizer.zero_grad()
        for batchNum in range(numtrainingBatches):
            b_labels, b_data = getMiniBatchData(y, M, train_minibatchIndices[batchNum], PRS_betas,  covars_all)
            # convert data to torch & move it to CUDA if available
            b_data = dataToDevice(b_data, device)
            b_labels = torch.from_numpy(b_labels).to(device)
            b_labels = b_labels.view(b_data.shape[0],-1 )  # make it the same shape as output, use b_labels instead of b_labels.shape[0], to ensure we have the correct 2nd dimension!

           # yhat.shape
           # b_labels.shape
           # b_labels_bad.shape
           # yhat_bad.shape
           # yhat.shape
            #loss_bad

     


            # Forward Propagate
            yhat = model(b_data)
            
            # perform full learning cycle, FP, BP and update weights
            #optimizer.zero_grad()   # zero the gradient buffers

            # loss function
            loss = criterion(yhat, b_labels)   / accumulation_steps

            # addRegularizerCosts(model, loss, args) # apply regularization to loss ( DONT add this, as we have this already in the optimizzer level weight_decay= l2)

            if t > 0 : # first epoch will want to evaluate loss based on no training
                loss.backward()   # Backprop
                
                # Gradient accumulation: only update params every accumulation_steps
                if (currentBatchNum+1) % accumulation_steps == 0:   
                    #torch.nn.utils.clip_grad_norm_(getModel(model).parameters(), 2) #  0.001
                    optimizer.step()                            # Now we can do an optimizer step
                    optimizer.zero_grad()                           # Reset gradients tensors

            epoch_trainingLoss+= loss.item() * (b_labels.shape[0] * b_labels.shape[1]) * accumulation_steps # get total loss # batch size: the criterion returns the AVG loss for the whole epoch, if we want to total loss, we need to multiply by batch_size: https://discuss.pytorch.org/t/mnist-dataset-why-is-my-loss-so-high-beginner/62670

            # cosmetics: update prograss bar
            barPos =  float(currentBatchNum) / numtrainingBatches # the bar position in %
            barPos = round(20 * barPos) # scale it to 1-20
            if suppressPrint == False : 
                sys.stdout.write('\r')
                sys.stdout.write("[%-20s] %d%%" % ('='*barPos, 5*barPos))
                sys.stdout.flush()  
            currentBatchNum += 1
            ###################### end of Training ########################


        # Evaluate Training (we now have to do this every turn, as we save every epoch's weights so that we dont have to retrain)
        epoch_trainingLoss= epoch_trainingLoss /numtrainingBatches # divide by number of batches to get a loss comparable between training/valid
        out_str =  out_str + " / Train loss: "+ str( round( epoch_trainingLoss, 3) )   #   + " / R^2: " + str(  round( train_accuracy,4) ) + " / S: " + str(  round( train_spearmanr[0],4) ) + "(" + str(  round( train_spearmanr[1],4) )  +")"
        results["train_loss"].append(epoch_trainingLoss)  
        

        # Eval Test 
        torch.set_grad_enabled(False) ; model.eval() # disable gradients to save memory
        epoch_validLoss = 0.
        numTestBatches =    len(test_minibatchIndices)
        for batchNum in range(numTestBatches):
        #for b_data, b_labels in zip(test_X, test_Y):
            b_labels, b_data = getMiniBatchData(y, M, test_minibatchIndices[batchNum], PRS_betas,  covars_all)
            b_data = dataToDevice(b_data, device) 
            b_labels = torch.from_numpy(b_labels).to(device) 
            b_labels = b_labels.view(b_data.shape[0],-1 )  # make it the same shape as output, use b_labels instead of b_labels.shape[0], to ensure we have the correct 2nd dimension!
       
            # Forward Propagate, based on what data was supplied
            yhat = model(b_data) # need to compute the Yhat again, as this is the yhat AFTER updating the weights, not before as in 'learn()' function

            #yhat = model(b_data) 
            loss = criterion(yhat, b_labels)

            # get total loss for iteration
            epoch_validLoss+= loss.item() * (b_labels.shape[0] * b_labels.shape[1]) # batch size: the criterion returns the AVG loss for the whole epoch, if we want to total loss, we need to multiply by batch_size: https://discuss.pytorch.org/t/mnist-dataset-why-is-my-loss-so-high-beginner/62670
            b_data = None

        # Evaluate Validation
        # valid_accuracy, valid_spearmanr = getCoefficientOfDetermination(yhats_saved,b_labels_saved)
        epoch_validLoss= epoch_validLoss /numTestBatches # divide by number of batches to get a loss comparable between training/valid
        results["valid_loss"].append(epoch_validLoss)
        out_str =  out_str + " | Valid loss: "+  str(round(epoch_validLoss,3))# + " / R^2: " +  str(round( valid_accuracy,4)) + " / S: " +  str(  round( valid_spearmanr[0],4) )  + "(" +  str(  round( valid_spearmanr[1],4) ) + ")"
       
        # if training has improved over the best so far, reset counter
        if results['lowestLoss'] is None or epoch_validLoss < results['lowestLoss']:
            results['lowestLoss'] = epoch_validLoss
            results['lowestLoss_epoch'] = t
            training_hasntmprovedNumIts = 0 
        else :
            training_hasntmprovedNumIts += 1

        # save the entire model weights for this epoch
        if saveModelLocation is not None:
            torch.save(getModel(model).state_dict(), saveModelLocation + str(t))  # could alternatively store it on the CPU as:  for k, v in state_dict.items():  state_dict[k] = v.cpu()   # https://discuss.pytorch.org/t/how-to-get-a-cpu-state-dict/24712
            torch.save(optimizer.state_dict(), saveModelLocation + str(t) + "_optimizer")  # also save the optimizer, otherwise we won't be able to resume training with a 100% accuracy, as optimizer has momentum etc which affect LR    # https://discuss.pytorch.org/t/loading-optimizer-dict-starts-training-from-initial-lr/36328      # and  https://discuss.pytorch.org/t/discontinuity-in-learning-rate-value-when-resuming-training-from-checkpoint/93128   
            writeOutEpochStats(saveModelLocation, results, t, training_hasntmprovedNumIts)
    
        elapsed_time = time.time() - start_time 
        if suppressPrint == False : print(out_str + " / " + str( round(elapsed_time) ) + " secs (LR: " + str(optimizer.state_dict()['param_groups'][0]['lr'])  + ")" , flush=True)  # lr_then = optimizer.state_dict()['param_groups'][0]['lr'] # str(round(currentEta,5))
        
        # update learning rate
        if decayRate != -1 and t > 0 :
            lr_scheduler.step()

        
        t += 1

    setModelMode(model, model_training_orig)
    return ( { "results" : results})  


###############################################################################
# Utils:
#region Utils:
###############################################################################
def writeOutEpochStats(saveModelLocation, results, t, training_hasntmprovedNumIts) :

    # write out the epoch itself and the number of epochs since training hasn't improved
    with open(saveModelLocation + str(t) +"_stats", "w") as file: 
        file.write(str(t) + "\t" +str(training_hasntmprovedNumIts) )
        
    with open(saveModelLocation + str(t) + "_results", "w") as file: 
        file.write( str(results['lowestLoss']) + "\n")      
        file.write( str(results['lowestLoss_epoch']) + "\n")  
        
        line=""
        if(len(results["epochs"]) > 0) :
            line=str(results["epochs"][0])
            for i in range(1,len(results["epochs"])) :
                line+="," + str(results["epochs"][i])
        file.write(line + "\n" ) 
        
        line=""
        if(len(results["train_loss"]) > 0) :
            line=str(results["train_loss"][0])
            for i in range(1,len(results["train_loss"])) :
                line+="," + str(results["train_loss"][i])
        file.write(line + "\n" ) 
        
        line=""
        if(len(results["valid_loss"]) > 0) :
            line=str(results["valid_loss"][0])
            for i in range(1,len(results["valid_loss"])) :
                line+="," + str(results["valid_loss"][i])
        file.write(line + "\n" ) 

    # results["epochs"] = [0,1,2,3]
    # results["train_loss"]  = [0.5,1.6,2.9,0.3]
    # results["valid_loss"]  = [0.5,1.6,2.9,0.3]
    # results['lowestLoss'] = 2.9
    # results['lowestLoss_epoch'] = 2
    
    # results2 = {}
    # results2["epochs"] = [0,1,2,3]
    # results2["train_loss"]  = [0.5,1.6,2.9,0.3]
    # results2["valid_loss"]  = [0.5,1.6,2.9,0.3]
    
    # results2['lowestLoss'] = 2.9
    # results2['lowestLoss_epoch'] = 2

#location = saveModelLocation+ "_" + str(t)
def loadEpochStats(location) :
    results = {}
    results["epochs"] = list()
    results["train_loss"]  = list()
    results["valid_loss"]  = list()
    results['lowestLoss'] = None
    results['lowestLoss_epoch'] = -1

    with open(location +"_stats", "r") as statsFile:
        itmp = statsFile.readline().rstrip().split()
        t = int(itmp[0])
        training_hasntmprovedNumIts = int(itmp[1])

    

    with open(location +"_results", "r") as resultsFile:
            counter = 0
            for i in resultsFile:
                if counter == 0 : 
                    x = i.rstrip()
                    x = None if x == 'None' else float(x)
                    results['lowestLoss'] = x
                    #print(results['lowestLoss'] )
                elif counter == 1 : 
                    results['lowestLoss_epoch'] = int(i.rstrip())
                    #print(results['lowestLoss_epoch'] )
                else :
                    itmp = i.rstrip().split(",")
                    if itmp[0] == '' : itmp = []

                    if counter == 2 :
                        itmp = [ int(x) for x in itmp ]
                        results["epochs"] = itmp
                        
                    elif counter == 3 :
                        itmp = [ float(x) for x in itmp ]
                        results["train_loss"] = itmp
                        
                    elif counter == 4 :
                        itmp = [ float(x) for x in itmp ]
                        results["valid_loss"] = itmp
                        
                counter = counter+1
        
    return results, t, training_hasntmprovedNumIts



def writeDebugOutput(debugOut,yhats_saved, b_labels_saved,b_weights_saved,t) :
#    yhats_all = torch.cat(yhats_saved)
#    b_labels_saved = torch.cat(b_labels_saved)
#    yhats_all = yhats_all.cpu().numpy()
#    b_labels_saved = b_labels_saved.cpu().numpy()
#    

    
    if type(yhats_saved) is list :
        yhats_all = np.concatenate(yhats_saved)
        b_labels_all = np.concatenate(b_labels_saved)
        b_weights_all = np.concatenate(b_weights_saved)
    else :
        yhats_all = yhats_saved
        b_labels_all = b_labels_saved
        b_weights_all = b_weights_saved
    with open(debugOut+ "_" + str(t), "w") as file: 
        file.write("y(" + str(np.mean(b_labels_all)) +")"  + "\t" + "yhat" + "\t"+ "weight"  + "\n")
        for i in range(len(b_labels_all)) :
            file.write( str(b_labels_all[i])  + "\t" + str(yhats_all[i]) + "\t"+ str(b_weights_all[i])  + "\n")   
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

def getMiniBatchData(y, M, batchIndices, PRS_betas,  covars_all) : 
    # get SNP data
    if M is not None: # if there is SNP data
        b_data = M[batchIndices]

        # multiply out the SNPs with the PRS weights (if any) at float32 precision
        if PRS_betas is not None : 
            b_data = b_data * PRS_betas #   b_data @ PRS_betas: this would NOT be what we want, as this would create a vector of size (minibatch,1), ie the per indi PRS. But here we just want each SNP [0,1,2] to be multiplied by their respective GWAS weight. sum(b_data[0]) would be then recover the PRS  b_data2[0] ==  sum(b_data3[0])
        #else : 
        b_data = b_data.astype('float32') # still need to cast this as the above may create float64
    else : # if there is no SNP data, then we must have had covars data, so then the batch data will just be the covars
        b_data = covars_all[batchIndices]

    # depending on if we have covars, we concat that into the data (here b_data can ONLY come from M, as we just made sure above)
    if covars_all is not None and M is not None: 
        b_data_all = np.zeros( (b_data.shape[0],b_data.shape[1] + covars_all.shape[1] ), dtype = np.float32 ) # create blank 2D array that could accommodate both the SNP and covars data
        b_data_all[:,0:b_data.shape[1]] = b_data # paste in the SNP data in the first part
        b_data_all[:,b_data.shape[1]:b_data.shape[1] + covars_all.shape[1]] = covars_all[batchIndices] # paste in the covars data in the second half
        b_data = b_data_all
        del b_data_all; gc.collect();

    # get the phenotypes
    b_labels = y[batchIndices]


    return b_labels, b_data



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




def dataToDevice(b_data, device) : # cast data to the device, if this is a 'tree' type NN, then the data will be a list of arrays, so we cannot directly cast it to cude all at once, but have to loop through and do it one by one
    
    #b_data = b_data.astype(np.float32) # need to enforce this. No actually we do this in the getminibatchdata now
    b_data =  torch.from_numpy(b_data).to(device) 
 
    #if self.half :              
    #    y_data = y_data.half()      
     #   b_data = b_data.half()


    return(b_data)

#endregion