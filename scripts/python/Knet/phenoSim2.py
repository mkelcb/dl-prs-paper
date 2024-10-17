# -*- coding: utf-8 -*-
"""
Created on Sun May  5 18:43:39 2019

@author: mk23
"""

# Given a UKBB bgen variant stats file, it gives back a range that equals a 500SNP flanking window around 2 variants

import argparse
import numpy as np
from itertools import combinations
import os
import struct
import time
import statsmodels.api as sm
import random
    
#from com.application.utils.regression import multiReg, uniReg
#from com.application.utils.geno_qc import removeList, genoQC_all, standardise_Genotypes, getSizeInMBs 

from knet_IO import loadPLINK, loadMatrixFromDisk, writePLINK_Pheno, loadPLINKPheno





#34270 / (34270 + 137088)
def printElapsedTime(start,end, text ="") : # https://stackoverflow.com/questions/27779677/how-to-format-elapsed-time-from-seconds-to-hours-minutes-seconds-and-milliseco
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print(text + "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds), flush=True)
        


def get_d_code(dataType) :     # need to  convert dataType into single letter codes
    if( dataType == 'float32') : d_code = 'f'
    elif ( dataType == 'float16') : d_code = 'e' 
    else : d_code = 'd' 
    
    return(d_code)
    



# Sum over an axis is a reduction operation so the specified axis disappears.  
def zscore(a, axis=0):
    a = a.astype('float32')  # if we don't cast this then this would upconvert everything to float64
    mns = a.mean(axis=axis)
    sstd = a.std(axis=axis)
    a = (a - mns) / sstd
    
    return a, mns, sstd

def writeSimPhenoToDisk(outFile, y) :
    with open(outFile , "w") as file: 
        for i in range(0, len( y ) ):
            file.write(str(y[i,0]) + "\n")
            



###############################################################################
# Data Load
###############################################################################   
    
def dataLoad(dataLoc, loadPhenos = True, replaceMissing = False) : # generic input data loading function that works for both PLINK binaries as well as LDPred PRS binaries
    # decide on if its a PLINK file, based on the extension: .bed file
    if os.path.exists(dataLoc + ".bed") :
        print("KNET: loading plink file")
        genotypeData = loadPLINK(dataLoc, loadPhenos = loadPhenos, replaceMissing = replaceMissing) 
        genotypeData["plink"] = True
    else : # create a dictionary with the exam same signature as the ones from the PLINK one
        print("KNET: loading PRS binaries")
        genotypeData = {}
        genotypeData["plink"] = False
        
        id_list = list()
        id_list.append( list())
        id_list.append( list())  
        with open(dataLoc + "_indis", "r") as id: # and indis file is created along the binaries, that has two columns: IID IID
            for i in id:
                itmp = i.rstrip().split()
                id_list[0].append( itmp[0] )
                id_list[1].append( itmp[1] )  
        genotypeData["IDs"] = id_list
        
        #genotypeData = loadPLINK("C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/data/genetic/short")
        genotypeData["M"] = loadMatrixFromDisk(dataLoc) 
        genotypeData["rsid"] = loadRegionsData(dataLoc +"_regions") # instead of SNP RSIds we will have regions like chrom_start_end, but these will correspond to the columns of the PRS matrix the same way as SNP IDs
  
    return(genotypeData)


#def dataLoad(dataLoc, loadPhenos = True, replaceMissing = False) : # generic input data loading function that works for both PLINK binaries as well as LDPred PRS binaries
#    # decide on if its a PLINK file, based on the extension: .bed file
#    if os.path.exists(dataLoc + ".bed") :
#        print("KNET: loading plink file")
#        genotypeData = loadPLINK(dataLoc, loadPhenos = loadPhenos, replaceMissing = replaceMissing) 
#        genotypeData["plink"] = True
#    else : # create a dictionary with the exam same signature as the ones from the PLINK one
#        print("KNET: loading PRS binaries")
#        genotypeData = {}
#        genotypeData["plink"] = False
#        
#        id_list = list()
#        id_list.append( list())
#        id_list.append( list())  
#        with open(dataLoc + "_indis", "r") as id: # and indis file is created along the binaries, that has two columns: IID IID
#            for i in id:
#                itmp = i.rstrip().split()
#                id_list[0].append( itmp[0] )
#                id_list[1].append( itmp[1] )  
#        genotypeData["IDs"] = id_list
#        
#        #genotypeData = loadPLINK("C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/data/genetic/short")
#        genotypeData["M"] = loadMatrixFromDisk(dataLoc) 
#        genotypeData["rsid"] = loadRegionsData(dataLoc +"_regions") # instead of SNP RSIds we will have regions like chrom_start_end, but these will correspond to the columns of the PRS matrix the same way as SNP IDs
#  
#    return(genotypeData)


###############################################################################
# Sim
###############################################################################



# lowmem version, we just pick interactions randomly from a pool until we haven't got enough (this is a 100x faster than the regular version)


#p = 50
#percCausal = 0.95
#additive = True
#orderInteraction = 2
#numInteractions = -1
#minOrderInteraction =2
#seed = 42
#h2 = 0.05
#
#genotypeData = loadPLINK("C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/data/genetic/short", loadPhenos = False, replaceMissing = True) 
#X = genotypeData["M"]
#
#X = X[0:2000,0:p]
###########################################################

#orderInteraction = minOrderInteraction = 4
# generates a genetic architecture based on various parameters
def buildGeneticArchitecture(p, percCausal = 1., additive = True, orderInteraction = 2, numInteractions = -1, minOrderInteraction =2, seed = 42) :
    np.random.seed(seed)
    random.seed(seed)
    
    # find out how many causals we will pick from that
    numCausals = int( np.floor(p * percCausal) )
    
    # get that many indices (save these)
    indices = np.asarray( range(p) )
    random.shuffle(indices)
    indices = indices[0:numCausals] # pick as many indices as we needed

    Betas_additive = np.array([])
   

    # if additive effects are enabled, we assign an effect size for all of these
    if additive :
        Betas_additive = np.random.normal( size=len(indices) ).astype(np.float32)
    

       
    #print("numCausals is: "   + str(numCausals))
    # if orderInteraction is not -1, we pick percCausal interactions of order orderInteraction: IE if we have to have 500 causals, with order 3, then we pick 500 times 3 indices (should I remove interaction with itself??)
    interactions = list()
    Betas_interaction_all = list()
    interactingPredictors = list()
    Betas_interaction_all_arrays = list() # conver the list of list to a list of arrays
    if orderInteraction != -1 :
        if numInteractions == -1 : numInteractions = numCausals # if we didn't specify the number of interactions we just use the number of causals as a basis
   

        for k in range(minOrderInteraction,orderInteraction+1) : # go through each level of interaction
            print("at order " + str(k) + " looking for numInteractions: " + str(numInteractions))
            start = time.time()
            
            currenOrderInteractions = list()
            Betas_interaction = list()

            failedTrials = 0
            #z = 0
            
            len_currenOrderInteractions = len(currenOrderInteractions)
            while len_currenOrderInteractions < numInteractions and failedTrials < 100: # keep going until we dont have enough interactions
                
                for i in range(len(currenOrderInteractions) , numInteractions) :
                    indices_working = indices.copy() # duplicate this
                    #print("attempt " + str(z)) ; z += 1
                    random.shuffle(indices_working) # randomise the remainder  
                    interactionIndices = indices_working[0:k] # pick the sufficient number for current level
                    interactionIndices.sort()
                    
                    currenOrderInteractions.append( interactionIndices )
                    
                    
                tempArray = np.array(currenOrderInteractions)
                unique_rows = np.unique(tempArray, axis=0)
                if (unique_rows.shape[0] < len(currenOrderInteractions)) : 
                    print("found duplicates: " + str(len(currenOrderInteractions) - unique_rows.shape[0]))
                    failedTrials += 1
                    
                currenOrderInteractions = list(unique_rows) # convert back to a list of arrays
                len_currenOrderInteractions = len(currenOrderInteractions)
                
            interactions.append(currenOrderInteractions) # have to add it here, as we are reassigning the variable above
            
            for i in range(0, len(currenOrderInteractions)) :
                interactionIndices = currenOrderInteractions[i]

                
                Beta =  np.float32( np.random.normal() ) # as this is one predictor actually just use one Beta obviously...
                Betas_interaction.append(Beta)


            Betas_interaction_all.append(Betas_interaction)
            #numInteractions = numInteractions*3 # we always triple the number of interactions
            timeTook =  round(time.time() - start)
            print("took length: " + str(timeTook) )
      

        for i in range(0,len(Betas_interaction_all)) :
            Betas_interaction = np.asarray(Betas_interaction_all[i], dtype=np.float32)
            Betas_interaction_all_arrays.append(Betas_interaction)
            
 
        # find the list of predictors that were involved in interactions, this may NOT be the full list if ''causals' we picked
        interactingPredictors = np.array([], dtype=np.int)
        for i in range(0,len(interactions)) :
            newElements = np.concatenate(np.array( interactions[i] ))
            interactingPredictors =  np.concatenate( [interactingPredictors,newElements] )
            
        interactingPredictors = np.unique(  interactingPredictors )


    # if we didn't request an additive effect, then do NOT save the indices that were used
    if additive == False: indices = np.asarray([]).astype(np.int)
 
    
    # return dict of:  list of causal SNPs, and list of interactions
    sim = {'indices_additive':indices, 'Betas_additive':Betas_additive, 'interactions':interactions, 'interactingPredictors':interactingPredictors, 'Betas_interaction_all_arrays': Betas_interaction_all_arrays }
    return(sim )


# given a geneticArchitecture dic, it generates phenotypes
# The formula for simulating pheno with given narrow sense heritability for an additive phenotype,
# y = g + noise
# noise ~N(0, sqrt(1-h2) )  # the noise is all the variance that is NOT due to h2
# g = GV * scale # scale the pure genetic values 
# GV = sum(X * Beta) # the breeding values are simply the linear combination of the alleles and the effects
# scale = sqrt( h2 / var(GV) ) # scale is the ratio of the desired h2, and the variance of the pure genetic values
# it is not just var(Beta), because the final breeding values are a product of the genotypes too

# to simulate phenotype due to interactions, with the appearance of given narrow sense h2:
# y' = g' + noise
# g' = GV * scale' # scale the pure genetic values with a scale proportionate to the apparent additive effect of variants
# scale = sqrt(h2 / var(GV') ) # variance of altered genetic values, sqrt to bring it back to the same scale as the phenotype
# GV = sum(X * Beta) 
# GV' <- lm(GV ~ X) 
# - that is, I fit a multi OLS to obtain additive effects that could have produced 'GV', and then take the fitted values from that model


def generateSim_AdditivePhenoWith_h2(X, sim, h2, seed = 42) : # generates additive phenotypes
    print ("SIMULATING TRUE ADDITIVE PHENOS")
    np.random.seed(seed)
    random.seed(seed)
    n = X.shape[0]

    indices = sim['indices_additive']
    Betas_additive = sim['Betas_additive']

    # simulate noise with the correct SD (1-h2), IE everything that is not heritable is due to noise
    noise = np.random.normal(size=n, scale=np.sqrt(1.0 - h2) ).astype(np.float32)

    # create 'pure' genetic value effects,
    GV = np.zeros( n, dtype =np.float32)   
    for i in range(0,len(indices)) :
        GV += (X[:, indices[i] ] * Betas_additive[i]  )
    
    var_gv = np.var(GV)
    print ("h2:", h2," / var_gv:", var_gv )
    scale = np.sqrt( h2 /var_gv ) # get the scale
    g = GV * scale # scale the pure genetic values to be proportionate to what is expected from h2


    # the final phenotype is a sum of the true genetic variance and the noise
    y = g + noise
    
    # normalize 
    y = zscore(y)[0]
    
    # reshape to correct dimensions
    y = y.reshape(len(y),1)

    return(y)

# simulate a phenotype with a given level of H2 from interactions (this is H2, not h2, as the total genetic variance will include the interaction effects)
def generateSim_AdditivePhenoWith_h2_naive(X, sim, h2, seed = 42) : # generates additive phenotypes
    print("simulating pheno with broad sense H2")
    np.random.seed(seed)
    random.seed(seed)
    n = X.shape[0]

    indices = sim['indices_additive']
    interactions = sim['interactions']
    Betas_additive = sim['Betas_additive']
    interactingPredictors = sim['interactingPredictors']
    Betas_interaction_all_arrays = sim['Betas_interaction_all_arrays']

    # obtain a total list of SNPs that are involved in interactions OR additive effects
    allCausalPredictors =  np.concatenate( [indices,interactingPredictors] )
    allCausalPredictors = np.unique(  allCausalPredictors ).astype(np.int32)
        
    # simulate noise with the correct SD (1-h2), IE everything that is not heritable is due to noise
    noise = np.random.normal(size=n, scale=np.sqrt(1.0 - h2) ).astype(np.float32)

    # find out the minimum and maximum order of interactions
    orderInteraction = -1
    if len(interactions) != 0 : 
        orderInteraction = len(interactions[-1][0]) # the max order interaction, find last interaction and look at its length

    # create 'pure' genetic value effects,
   
    # Additive part
    GV = np.zeros( n, dtype =np.float32)   
    for i in range(0,len(indices)) : # add additives
        GV += (X[:, indices[i] ] * Betas_additive[i]  )

    # interaction part
    if orderInteraction != -1 :
        for k in range(len(interactions)) : # go through each level of interaction
            currenOrderInteractions = interactions[k]
            Betas_interaction = Betas_interaction_all_arrays[k]

            for i in range(len(currenOrderInteractions)) :
                interactionIndices = currenOrderInteractions[i]
                Beta = Betas_interaction[i]  # as this is one predictor actually just use one Beta obviously...
                products = np.product(X[:,interactionIndices ], axis = 1)
                GV += ( products * Beta)

    
    var_gv = np.var(GV)
    print ("h2:", h2," / var_gv:", var_gv )
    
    
    scale = np.sqrt( h2 /var_gv ) # get the scale # get the scale, but this time use the 'apparent' genetic values
    g = GV * scale # scale the pure genetic values to be proportionate to what is expected from h2


    # the final phenotype is a sum of the true genetic variance and the noise
    y = g + noise
    
    # normalize 
    y = zscore(y)[0]
    
    # reshape to correct dimensions
    y = y.reshape(len(y),1)

    return(y)

    
    

def generateSim_AdditivePhenoWith_h2_interactions(X, sim, h2, seed = 42) : # generates a phenotype with a given apparent h2, but which is due to interactions
    np.random.seed(seed)
    random.seed(seed)
    n = X.shape[0]

    indices = sim['indices_additive']
    interactions = sim['interactions']
    Betas_additive = sim['Betas_additive']
    interactingPredictors = sim['interactingPredictors']
    Betas_interaction_all_arrays = sim['Betas_interaction_all_arrays']

    # obtain a total list of SNPs that are involved in interactions OR additive effects
    allCausalPredictors =  np.concatenate( [indices,interactingPredictors] )
    allCausalPredictors = np.unique(  allCausalPredictors ).astype(np.int32)
        
    # simulate noise with the correct SD (1-h2), IE everything that is not heritable is due to noise
    noise = np.random.normal(size=n, scale=np.sqrt(1.0 - h2) ).astype(np.float32)

    # find out the minimum and maximum order of interactions
    orderInteraction = -1
    if len(interactions) != 0 : 
        orderInteraction = len(interactions[-1][0]) # the max order interaction, find last interaction and look at its length

    # create 'pure' genetic value effects,
   
    # Additive part
    GV = np.zeros( n, dtype =np.float32)   
    for i in range(0,len(indices)) : # add additives
        GV += (X[:, indices[i] ] * Betas_additive[i]  )

    # interaction part
    if orderInteraction != -1 :
        for k in range(len(interactions)) : # go through each level of interaction
            currenOrderInteractions = interactions[k]
            Betas_interaction = Betas_interaction_all_arrays[k]

            for i in range(len(currenOrderInteractions)) :
                interactionIndices = currenOrderInteractions[i]
                Beta = Betas_interaction[i]  # as this is one predictor actually just use one Beta obviously...
                products = np.product(X[:,interactionIndices ], axis = 1)
                GV += ( products * Beta)

#####################################################################
    # Obtain per predictor estimated Betas by fitting linear model against the GV, these Betas try to approximate what the apparent additive effect of any SNP would be
    
#    # fitting via series of univariate regressions ( results in a poor fit, probably due to correlated data)
#    beta_additive_est = np.zeros( len(allCausalPredictors), dtype =np.float32)  
#    GV_prime = np.zeros( n, dtype =np.float32)  
#    for i in range(len(allCausalPredictors)) :
#        predictor_i = allCausalPredictors[i]
#        X_i = X[:,predictor_i]
#        X_withIntercept = sm.add_constant(X_i)
#        model = sm.OLS(GV, X_withIntercept)
#        results = model.fit()
#        beta_additive_est[i] = results.params[-1]
#        GV_prime += X[:,predictor_i] * beta_additive_est[i]
#  
#    fitQuality = np.corrcoef(GV_prime,GV)[0,1]    # additive = 0.68  oi2 = 0.67 oi3= 0.65, oi4: 0.47
#    np.corrcoef(beta_additive_est,Betas_additive)
#    beta_additive_est[2]
#    Betas_additive[2]
    
        
    # Multiple regression,
    X_ = X[:,allCausalPredictors] # extract out from the original full predictor matrix the ones that are involved in causal effects
    X_ = sm.add_constant(X_) # add intercept
    model = sm.OLS(GV, X_) # fit multireg model against the true genetic values
    results = model.fit()
    beta_additive_est = results.params[1:len(results.params)] # dont want the intercept
    GV_prime = results.fittedvalues # model fit values GV'
        
    # test Fit: it is pretty tight for additive ~1, but detoriates for higher order interactions
    fitQuality = np.corrcoef(GV_prime,GV)[0,1]   # additive: 1, oi2: 0.8, oi3= 0.7, oi4: 0.52 
    print("additive estimation fit quality for order interaction",orderInteraction ,":", fitQuality)
#    np.corrcoef(beta_additive_est,Betas_additive)
#    beta_additive_est[3]
#    Betas_additive[3]
#################################################################
      
    
    var_gv = np.var(GV_prime)
    print ("h2:", h2," / var_gv:", var_gv )
    
    
    scale = np.sqrt( h2 /var_gv ) # get the scale # get the scale, but this time use the 'apparent' genetic values
    g = GV * scale # scale the pure genetic values to be proportionate to what is expected from h2


    # the final phenotype is a sum of the true genetic variance and the noise
    y = g + noise
    
    # normalize 
    y = zscore(y)[0]
    
    # reshape to correct dimensions
    y = y.reshape(len(y),1)

    return(y)











def writeOutSimDetails(outFile, sim, minOrderInteraction = 2) :

  #  writePLINK_Pheno(outFile,  sim['y'])

    with open(outFile + ".indices", "w") as file: 
        for i in range(0, len( sim['indices_additive'] ) ):
            file.write(str(sim['indices_additive'][i]) + "\n")
            
            len( sim['Betas_additive'] )
            sim['Betas_additive']
            
    #if sim['Betas_additive'] is not None :
    with open(outFile + ".Betas_additive", "w") as file: 
        for i in range(0, len( sim['Betas_additive'] ) ):
            file.write(str(sim['Betas_additive'][i]) + "\n")

    with open(outFile + ".interactingPredictors", "w") as file: 
        for i in range(0, len( sim['interactingPredictors'] ) ):
            file.write(str(sim['interactingPredictors'][i]) + "\n")
            
    for j in range(0, len( sim['Betas_interaction_all_arrays'] ) ):     # nested structure, for each level of interaction there is a list of indices
        with open(outFile + ".Betas_interaction_oi_" + str(j+minOrderInteraction), "w") as file:  # +2 as there is a minimum
            for i in range(0, len( sim['Betas_interaction_all_arrays'][j] ) ):
                file.write(str(sim['Betas_interaction_all_arrays'][j][i]) + "\n")
    
    for j in range(0, len( sim['interactions'] ) ):     # nested structure, for each level of interaction there is a list of indices
        with open(outFile + ".interactions_oi_" + str(j+minOrderInteraction), "w") as file:  # +minOrderInteraction as there is a minimum
            for i in range(0, len( sim['interactions'][j] ) ):
                interactingIndices = str(sim['interactions'][j][i][0]) # add the first index  for formatting purposes
                
                for k in range(1, len(sim['interactions'][j][i])):
                    interactingIndices = interactingIndices + "\t" + str(sim['interactions'][j][i][k])
                    
                file.write(interactingIndices + "\n")
    
    
    
    
def readSimDetails(outFile, orderInteraction=4, minOrderInteraction = 2) :
 
   # y = loadPLINKPheno(outFile + ".pheno", caseControl = False)
   # y = y[:,np.newaxis]
    
    
    indices = list()
    with open(outFile + ".indices", "r") as file:
        for i in file:
            itmp = i.rstrip().split()
            indices.append( int(itmp[0]) )
    indices = np.array(indices)


    Betas_additive = list()
    with open(outFile + ".Betas_additive", "r") as file:
        for i in file:
            itmp = i.rstrip().split()
            Betas_additive.append( np.float32(itmp[0]) )
    Betas_additive = np.array(Betas_additive)



    interactingPredictors = list()
    with open(outFile + ".interactingPredictors", "r") as file:
        for i in file:
            itmp = i.rstrip().split()
            interactingPredictors.append( int(itmp[0]) )
    interactingPredictors = np.array(interactingPredictors)
    


    Betas_interaction_all_arrays = list()
    for j in range(minOrderInteraction, orderInteraction+1 ):     # nested structure, for each level of interaction there is a list of indices

        if os.path.exists(outFile + ".Betas_interaction_oi_" + str(j)) :
            with open(outFile + ".Betas_interaction_oi_" + str(j), "r") as file:  # +2 as there is a minimum
                Betas_interaction = list()
                for i in file:
                    itmp = i.rstrip().split()
                    Betas_interaction.append(  np.float32(itmp[0]) )
                    
                Betas_interaction = np.array(Betas_interaction)
                Betas_interaction_all_arrays.append(Betas_interaction)
            
            
            
    interactions = list()
    for j in range(minOrderInteraction, orderInteraction+1 ):     # nested structure, for each level of interaction there is a list of indices

        if os.path.exists(outFile + ".interactions_oi_" + str(j)) :
            with open(outFile + ".interactions_oi_" + str(j), "r") as file:  # +2 as there is a minimum
                interactions_order = list()
                for i in file:
                    itmp = i.rstrip().split()
                    interactingIndices = list()
                    for z in range(0, len(itmp)) :
                        
                        interactingIndices.append(  int(itmp[z]) )
                    
                    interactingIndices = tuple(np.array(interactingIndices) )
                    interactions_order.append(interactingIndices)
     
                
                interactions.append(interactions_order)
            
    
    sim_loaded = {'indices_additive':indices, 'interactions':interactions, 'interactingPredictors':interactingPredictors, 'Betas_additive':Betas_additive, 'Betas_interaction_all_arrays': Betas_interaction_all_arrays }
    
    return(sim_loaded )

###############################################################################


def loadRegionsData(regions) :
    # loads in regions data file line by line
    regionsData = list()
    with open(regions, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            regionsData.append( itmp[0] )
    
    return(regionsData)


def loadMatrixFromDisk(location, dataType ="float32") :
    d_code = get_d_code(dataType)
    
    # load id file to get dimensions
    with open(location + ".id", "r") as idFile:
        itmp = idFile.readline().rstrip().split()
        nrows = int(itmp[0])
        ncols = int(itmp[1])
        
    # how many elements to expect in the binary in total
    totalNum =nrows * ncols
    
    # open binary file
    with open(location + ".bin", "rb") as BinFile:
        BinFileContent = BinFile.read()
  
    # reformat data into correct dimensions
    flat = np.array( struct.unpack(d_code*totalNum, BinFileContent  ), dtype = dataType )
    data = flat.reshape(nrows,ncols)
    return(data)

    


##################################################################
# h2=0.05
# p = numVars
# out = 'C:/softwares/Cluster/GIANT/phenoSim/sim' 
# plinkFile1 = "C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/data/genetic/short"
def generatePhens(out,plinkFile1, numOrders, numVersions = 10,plinkFile2 = None,plinkFile3 = None, percCausal = 1., seed = 42, h2 = 0.05, numInteractions = -1) :

    # find out how many predictors the genotype has: just check the number of lines in either the .bim or the _regions file
    if os.path.exists(plinkFile1 + ".bim") : variantsFile = plinkFile1  + ".bim"
    else : variantsFile = plinkFile1  + "_regions"

    numVars = 0
    with open(variantsFile, "r") as id:
        for i in id:
            numVars += 1

    sims = list()
    for k in range(1, numOrders +1) :
        print("simulating phenos for order", k)
        
        sims_versions = list()
        for c in range(numVersions) : # also make 10 different versions of the same pheno
            if k == 1 :
                print("simulating additive pheno")
                sim = buildGeneticArchitecture(numVars, percCausal = percCausal, additive = True, orderInteraction = -1, numInteractions = numInteractions, minOrderInteraction = -1, seed = seed+c)
              
            else :
                sim = buildGeneticArchitecture(numVars, percCausal = percCausal, additive = False, orderInteraction = k, numInteractions = numInteractions, minOrderInteraction = k, seed = seed+c) 
            sims_versions.append(sim)
            
            # write sim to disk
            writeOutSimDetails(out  +"_" +str(k) + "_v" + str(c), sim, minOrderInteraction = k)
        sims.append(  sims_versions  )  

    # write out all 10 versions of the phenotype all at once
    outPutPhenos(out+"_train",plinkFile1,sims,h2, seed = seed +1 )
    if plinkFile2 != None : outPutPhenos(out +"_valid",plinkFile2,sims,h2, seed = seed + 2)
    if plinkFile3 != None : outPutPhenos(out +"_test",plinkFile3,sims,h2, seed = seed + 3)
    
        
 
def outPutPhenos(out,plinkFile1, sims, h2, seed = 42) :
    print("loading PLINK from is: " , plinkFile1)
   
    # this must be PLINK data, but it it CAN be in PRS format (IE a plink file coverted to PRS)
    genotypeData = dataLoad(plinkFile1, loadPhenos = False) # loadPLINK(args.knet, loadPhenos = False) 
    #plinkVersion = genotypeData["plink"]
    X = genotypeData["M"]
    #irsIds = genotypeData["rsid"]
    #IDs = genotypeData["IDs"] 
          
  
    for i in range(len(sims)) : # go throug the orders of interaction
        sims_versions = sims[i]
        
        for c in range(len(sims_versions)) : # go throug the versions
            sim = sims_versions[c]
            interactions = sim['interactions']
            orderInteraction = 1

            # find out the minimum and maximum order of interactions
            if len(interactions) != 0 : orderInteraction = len(interactions[-1][0]) # the max order interaction, find last interaction and look at its length
    
            
            if orderInteraction == 1: y = generateSim_AdditivePhenoWith_h2(X, sim, h2, seed = seed + i)
            else : y = generateSim_AdditivePhenoWith_h2_interactions(X, sim, h2, seed = seed + i) # we dont want the same noise for all  phenos, so we add i    
            
            outFile = out +"_" +str(orderInteraction) + "_v" + str(c)
            writeSimPhenoToDisk(outFile,y)
        # write pheno to disk
    
    
    print("simulations to:" , out )









#955*0.95 * 2

def runMain(args) :
    print('generating phenos from  ', args.plinkFile1, "numOrders:",args.numOrders , "numVersions:",args.numVersions, " with h2:", args.h2, " % causal:", args.percCausal, " num Interactions:", args.numInteractions)    
    out = args.out
    plinkFile1 = args.plinkFile1
    numOrders = args.numOrders
    numVersions = args.numVersions
    percCausal = args.percCausal
    
    h2 = args.h2
    
    seed = args.seed
    plinkFile2 = None
    plinkFile3 = None
    if  args.plinkFile2 != "0" : plinkFile2 = args.plinkFile2
    if  args.plinkFile3 != "0" : plinkFile3 = args.plinkFile3
    
    
    generatePhens(out,plinkFile1, numOrders, numVersions,plinkFile2, plinkFile3, percCausal, seed,h2,args.numInteractions )





if __name__ == '__main__':   
    parser = argparse.ArgumentParser()

    parser.add_argument("--out",required=True, help='an output location is always required')
    parser.add_argument("--plinkFile1",required=True, help='Location of the PLINK or PRS file.')  
    parser.add_argument("--plinkFile2",required=False, help='(optional) Location of a second PLINK or PRS file.', default="0")  
    parser.add_argument("--plinkFile3",required=False, help='(optional) Location of a third PLINK or PRS file.', default="0")  

    parser.add_argument("--numOrders",required=False, help='from 1 to X, the orders of interaction phenotypes to generate  (default=4)', default=4, type= int) 

    parser.add_argument("--numVersions",required=False, help='how many versions of each phenotype to generate (default=20)', default=20, type= int) 

    parser.add_argument("--percCausal",required=False, help='percentage of causal vars', default=0.5, type= float) 
    parser.add_argument("--seed",required=False, help='random seed', default=42, type= int) 

    parser.add_argument("--h2",required=False, help='narrow sense heritability', default=0.5, type= float) 

    parser.add_argument("--numInteractions",required=False, help='The number of interactions. -1 to set it equal to the number of causals', default=-1, type= int) 



    parser.set_defaults(func=runMain)
        
    args = parser.parse_args()
    args.func(args)



