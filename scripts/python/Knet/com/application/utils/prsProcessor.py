# -*- coding: utf-8 -*-
"""
Created on Sun May  5 18:43:39 2019

@author: mk23
"""

# Given a UKBB bgen variant stats file, it gives back a range that equals a 500SNP flanking window around 2 variants

import argparse
import numpy as np

from itertools import combinations

from regression import multiReg, uniReg
import struct


from pathlib import Path


def findMatchingIndices (fullList, sublist) :
      
    matching_indices = np.zeros(len(sublist), dtype=np.int32 )
    for i in range( len(sublist) ) : matching_indices[i] = fullList.index(sublist[i])
    return matching_indices


def get_d_code(dataType) :     # need to  convert dataType into single letter codes
    if( dataType == 'float32') : d_code = 'f'
    elif ( dataType == 'float16') : d_code = 'e' 
    else : d_code = 'd' 
    
    return(d_code)
    
# writes a matrix onto disk in a binary (2 files, 1 that stores the dimensions, the other the binary data)
def writeMatrixToDisk(location,data, dataType ="float32") :
    d_code = get_d_code(dataType)
    
    # get dimensions of matrix
    nrows = data.shape[0]
    ncols = data.shape[1]
    
    # write the dimensions onto disk
    with open(location + ".id", "w") as idFile: 
        idFile.write(str(nrows) + "\t" +str(ncols) )
        
    # flatten matrix
    flat = data.ravel()

    flatData = struct.pack(d_code*len(flat),*flat  )
    with open(location + ".bin", "wb") as flat_File: 
        flat_File.write(flatData) 
    

#def loadRegionsData(regions) :
#    # loads in regions data file line by line
#    regionsData = list()
#    with open(regions, "r") as id:
#        for i in id:
#            itmp = i.rstrip()
#            regionsData.append( itmp )
#    
#    return(regionsData)
#    
    
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
    
    
def recodeCaseControlQuantitative(status) : # transcode binary phenotype into a quantitative phenotype ( IE just turn the 1/2 into 0.0 and 1.0 )
    y = np.zeros( (len(status))  )
    for i in range(len(status)) :
        y[i] =(status[i] -1)  
    
    return(y)
    

#A1_alleles = list(bim.a1)
def recodeCaseControlOneHot(status) :
    # transcode phenotype into 1 hot array
    y = np.zeros( (len(status),2)  )
    for i in range(len(status)) :
        result = np.array([0,0])
        result[ int(status[i] -1) ] =1 # transocde phenotype into one-of-k format: in plink case =2, control =1, so we set the element to 1 to at the inde-1 position
        y[i] =result  
    
    return(y)
    
def loadPLINKPheno(location, caseControl = True, recodeCaseControl = True) :
    status = list()
    #print("FUCK YOU")
    with open(location, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            status.append( np.float32(itmp[2]) )
 
    if caseControl : # if trait is binary
        if recodeCaseControl : y = recodeCaseControlQuantitative(status) # if we want to treat it as a quantitative trait
        else : y = recodeCaseControlOneHot(status)# otherwise use the 'softmax' format                                      
    else : y = np.array(status)# otherwise just leave it alone
    
    return (  y.astype('float32')  )

##################################################################
# yLoc = 'C:/softwares/Cluster/GIANT/miniPRS/miniPRS.pheno'
# subsetRegions = "C:/softwares/Cluster/GIANT/miniPRS/subset_regions"
    
# yLoc = 'C:/softwares/Cluster/GIANT/miniPRS/miniPRS.pheno'
# subsetRegions = "C:/softwares/Cluster/GIANT/miniPRS/subset_regions"
# fullPRSLoc = 'C:/softwares/Cluster/GIANT/miniPRS/miniPRS'
# subsetindis = 'C:/softwares/Cluster/GIANT/miniPRS/subset_indis' 
# out = 'C:/softwares/Cluster/GIANT/miniPRS/miniPRS_subset' 

def processPRS(fullPRSLoc, out, yLoc = None, subsetRegions = None, performInteractionTest = False, multireg = False) :


    # load original Data
    id_list = list()
    id_list.append( list())
    id_list.append( list())  
    with open(fullPRSLoc + "_indis", "r") as id: # and indis file is created along the binaries, that has two columns: IID IID
        for i in id:
            itmp = i.rstrip().split()
            id_list[0].append( itmp[0] )
            id_list[1].append( itmp[1] )  

    
    #genotypeData = loadPLINK("C:/Users/mk23/GoogleDrive_phd/PHD/Project/Implementation/data/genetic/short")
    M = loadMatrixFromDisk(fullPRSLoc) 
    rsid = loadRegionsData(fullPRSLoc +"_regions") # instead of SNP RSIds we will have regions like chrom_start_end, but these will correspond to the columns of the PRS matrix the same way as SNP IDs
  



    if subsetRegions is None :
        y = loadPLINKPheno(yLoc, caseControl = False)
        
        if performInteractionTest :
            # attempt to fit interaction tests
            
            print("epistasis requested")
            interactionResults = epistasisTest(M,y)
            # write results to disk
            with open(out + "_interactionResults", "w") as resultsFile: 
                resultsFile.write("SNP1" + "\t" + "SNP2" + "\t" + "p"  + "\t" +"beta" + "\n" )
                for i in range(interactionResults.shape[0])   : 
                    resultsFile.write(str(rsid[ int(interactionResults[i,0]) ]) + "\t" + str( rsid[ int(interactionResults[i,1]) ]) + "\t" + str(interactionResults[i,2])  + "\t"  + str(interactionResults[i,3]) + "\n" )  
                    
        
        # first build a total PRS (by summing al cols) and fit a model and also write out the data:
#        prsSum= np.sum(M, axis=1)
#        prsSum= prsSum.reshape(len(prsSum), 1)
#        
#        results_sum = multiReg(y, prsSum, bonfCorrect = False, addIntercept = True) 
#        
#        # write results to disk
#        with open(out + "_sumResults", "w") as resultsFile: 
#            resultsFile.write("predictor" + "\t" + "p"  + "\t" +"beta" + "\n" )
#            for i in range(1, len(results_sum["p_values"])  ) : # -1 as we dont care about the intercept 
#                resultsFile.write(str(rsid[i-1]) + "\t" + str(results_sum["p_values"][i])  + "\t"  + str(results_sum["beta"][i]) + "\n" )              
#        # write data to disk so it can also be repeated in R
#        with open(out + "_sumDesign", "w") as resultsFile: 
#            resultsFile.write("y" + "\t" + "prs_all"  + "\n")
#            for i in range(1, len(y)  ) : # -1 as we dont care about the intercept 
#                resultsFile.write(str(y[i]) + "\t" + str(prsSum[i][0])  + "\n" )  
#     

        # Fit GWAS regression (IE univariate regression of each marker separately) 
        try:
            resultsUni = uniReg(y, M, bonfCorrect = False) 
    
            with open(out + "_uniRegResults", "w") as resultsFile: 
                resultsFile.write("predictor" + "\t" + "p"  + "\t" +"beta" + "\n" )
                for i in range(1, len(resultsUni["p_values"])  ) : # -1 as we dont care about the intercept 
                    resultsFile.write(str(rsid[i-1]) + "\t" + str(resultsUni["p_values"][i])  + "\t"  + str(resultsUni["beta"][i]) + "\n" )
        except:
            print("bad design matrix for unireg")
    
        
        
        # Attempt Multiple regression
        if multireg :
            try:
                resultsMulti = multiReg(y, M, bonfCorrect = False, addIntercept = True)
                varExplained = 1 - resultsMulti["sse"]
                
                # write results to disk
                with open(out + "multiRegResults", "w") as resultsFile: 
                    resultsFile.write("predictor(1-RSS:" + str(varExplained) + ")\t" + "p"  + "\t" +"beta" + "\n" )
                    for i in range(1, len(resultsMulti["p_values"])  ) : # -1 as we dont care about the intercept 
                        resultsFile.write(str(rsid[i-1]) + "\t" + str(resultsMulti["p_values"][i])  + "\t"  + str(resultsMulti["beta"][i]) + "\n" )
            except:
                print("bad design matrix for multireg")
            
        

        
        
    else :
        print("subsetting regions to " ,subsetRegions)
        rsid_subset = loadRegionsData(subsetRegions) # load the predictors to be used
        
        # find the indices of the subset region in the full PRS regions file
        #matchingIndices = np.array([rsid.index(i) for i in rsid_subset] )
        matchingIndices = findMatchingIndices(rsid,rsid_subset)# np.array([rsid.index(i) for i in rsid_subset] )
        print("subset regions to ", len(matchingIndices), " out of ", M.shape[1] )
        prsMatrix = M[:,matchingIndices]
        del M # free up RAM

        # write it to disk
        # inline code to make it more memory efficient
        dataType ="float32"
        d_code = get_d_code(dataType)
        
        # get dimensions of matrix
        nrows = prsMatrix.shape[0]
        ncols = prsMatrix.shape[1]
        
        # write the dimensions onto disk
        with open(out + ".id", "w") as idFile: 
            idFile.write(str(nrows) + "\t" +str(ncols) )
            
        # flatten matrix
        flat = prsMatrix.ravel()
        del prsMatrix # free up RAM
    
    

        totalFlatLen= len(flat) 
        numParts = 10
        tenP = totalFlatLen//numParts # split it into 10 parts
        startIndex = 0
        for i in range(numParts) :  
            mode="ab"
            endIndex = startIndex + tenP
            if i == (numParts-1) : endIndex = totalFlatLen # the last part should cover up until the end
            if i == 0: mode = "wb" # for the first time we want it to start a new file
            
            print("startIndex:" , startIndex, "endIndex is:" , endIndex, " out of totalFlatLen: ", totalFlatLen)
            flatData = struct.pack(d_code*len(flat[startIndex:endIndex]),*flat[startIndex:endIndex]  )
            with open(out + ".bin", mode ) as flat_File: 
                flat_File.write(flatData) 
            
            startIndex = endIndex


            
        print("written PRS binary matrix to:" , out )


# k =2
# X = M[0:9,:]
def epistasisTest(X, y, orderInteraction = 2) : # produces a design matrix that includes all higher order interactions
    # figure out the dimensions of the final matrix, and pre-allocate it, so we can just copy into it
    indices = np.asarray( range(X.shape[1]) )
    

    k = 2 # only test 2 way interactions
    allInteractions = np.array( list(  combinations(indices, k) ) )
    interaction_results = np.zeros( (allInteractions.shape[0], 4), dtype = np.float32)
    for i in range(allInteractions.shape[0]) : # go through each interaction term
        print("epistasis test:", i, "/", allInteractions.shape[0])
        # fit a model for each term: y = SNP1 + SNP2 + SNP1*SNP2
        SNP1_index = allInteractions[i,0]
        SNP2_index =allInteractions[i,1]
        SNP1 = X[:,SNP1_index]
        SNP2 = X[:,SNP2_index]
        interactionTerm = np.product(X[:,allInteractions[i,:] ], axis = 1)
        X_currentInteractionDesignMatrix = np.c_[SNP1,SNP2,interactionTerm ]
        try:
            resultsMulti = multiReg(y, X_currentInteractionDesignMatrix, bonfCorrect = False, addIntercept = True)
            # extract interaction term p_val and Beta
            pval = resultsMulti['p_values'][-1]
            beta = resultsMulti['beta'][-1]
            interaction_results[i] = np.array([SNP1_index,SNP2_index,pval,beta])
        except:
            print("bad design matrix for SNPS", SNP1_index, "and", SNP2_index)
            interaction_results[i] = np.array([SNP1_index,SNP2_index,-1.,-1.])  

    return(interaction_results)
    

def produceXInteractionOLS_staircase(X, orderInteraction = 2) : # produces a design matrix that includes all higher order interactions
    # figure out the dimensions of the final matrix, and pre-allocate it, so we can just copy into it
    indices = np.asarray( range(X.shape[1]) )
    
    interactionTerms = list()
    totalNumExtraPredictors = 0 # how many more predictors we need to add to the design matrix
    for k in range(2,orderInteraction+1) : # go through each level of interaction
        
        allInteractions = np.array( list(  combinations(indices, k) ) )
        totalNumExtraPredictors+= allInteractions.shape[0]
        interactionTerms.append(allInteractions)
        print("order " + str(k) + " adds number of predictors: " + str(allInteractions.shape[0]) )
    
    
    X_design = np.zeros( (X.shape[0], X.shape[1] + totalNumExtraPredictors ) , dtype =np.float32)
    X_design[:,0:X.shape[1]] = X # copy over the original
    counter = X.shape[1]
    for j in range(0, len(interactionTerms)) :
        allInteractions = interactionTerms[j]
        # original size, then loop through each order interaction (or do we have just 1 level at once??)
        # this is basically all combinations (IE order does not matter)

        for i in range(0,allInteractions.shape[0]) :
            X_design[:,counter] = np.product(X[:,allInteractions[i]], axis = 1) 
            counter += 1

    return(X_design)


def runMain(args) :
    print('Processing PRS from ',  args.fullprsLoc)    
    out = args.out
    fullPRSLoc = args.fullprsLoc
    yLoc= None
    subsetRegions= None
    performInteractionTest = False
    multireg = False
    if  args.subsetregions != "0" : subsetRegions = args.subsetregions
    if  args.yloc != "0" : yLoc = args.yloc
    if  args.epistasis != "0" : performInteractionTest = True
    if  args.multireg != "0" : multireg = True
    
    

    processPRS(fullPRSLoc, out, yLoc, subsetRegions, performInteractionTest,multireg)
  #  (675 - 628) /2

if __name__ == '__main__':   
    parser = argparse.ArgumentParser()

    parser.add_argument("--out",required=True, help='an output location is always required')
    parser.add_argument("--fullprsLoc",required=True, help='The fulle PRS binary file')  
    parser.add_argument("--subsetregions",required=False, help='(Optional) a file describing the regions to be used to subset the full PRS file. ', default="0")  
    parser.add_argument("--yloc",required=False, help='(Optional) a phenotype file. Only needed if you want to run regression analyses. ', default="0")  
    parser.add_argument("--epistasis",required=False, help='(Optional) if a 2-way epistatic interaction test should be performed ', default="0")  
    parser.add_argument("--multireg",required=False, help='(Optional) if multiple regression OLS should be attempted', default="0")  
       
    
    parser.set_defaults(func=runMain)
        
    args = parser.parse_args()
    args.func(args)



