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


from pyplink import PyPlink



import numpy as np
import struct
from pathlib import Path
import collections
import gc




##################################################################################################
# PLINK Binary / Phenotype Load
##################################################################################################
#location= "C:/Users/leo/Google Drive/PHD/0Rotations/!0Implementation/tests/knet/py/data/toy/wtccc2_hg19_toy"
#results = loadPLINK(location, recodeCaseControl = False)
#location = args.bfile
#males = fam.gender == 1

#for genotypes in bed.iter_geno_marker(bim.pos):
#    male_genotypes = genotypes[males.values]

#male_genotypes = genotypes[males.values]
#location = args.out
# np.single
def get_d_code(dataType) :     # need to  convert dataType into single letter codes: https://numpy.org/doc/stable/reference/arrays.dtypes.html
    if( dataType == 'float32') : d_code = 'f'
    elif ( dataType == 'float16') : d_code = 'e' 
    elif ( dataType == 'int32') : d_code = 'i' 
    
    else : d_code = 'd' 
    
    return(d_code)

def writePLINK_Pheno(location, y, indList = None) :
    y = y.ravel()
    
    with open(location + ".pheno", "w") as file: 
        for i in range(0, len(y) ):
            if indList is None : file.write("IND"+ str(int(i+1)) + "\tIND" + str(int(i+1)) + "\t" +str(y[i]) + "\n")
            else : file.write(str(indList[i]) + "\t" +str(indList[i]) + "\t" +str(y[i]) + "\n")


def writePLINK(location, M, chrNum= 1, rsIds = None) :
    # load plink file
    bed = PyPlink(location, mode='w', bed_format='SNP-major')
    for i in range(M.shape[1]):
        #print("writing SNP", i, "to file")
        bed.write_genotypes(M[:,i])
    
    with open(location + ".fam", "w") as file: 
        for i in range( M.shape[0] ):
            file.write("IND"+ str(int(i+1)) + "\tIND" + str(int(i+1)) + "\t0\t0\t0\t0" + "\n") 

    with open(location + ".bim", "w") as file: 
        for i in range( M.shape[1] ):
            SNPId="SNP" + str(int(i))
            if rsIds is not None : SNPId = rsIds[i]
            file.write( str(chrNum)+"\t"+ SNPId + "\t0\t" + str( int( (i+1) * 1000)) + "\tA\tB"  + "\n") 


#location = args.knet, loadPhenos = False
#male_genotypes = genotypes[males.values]
def loadPLINK(location, loadPhenos = True, caseControl = True, recodeCaseControl = True, replaceMissing = False) :
    # load plink file
    bed = PyPlink(location, mode='r', bed_format='SNP-major')
    bim = bed.get_bim()
    fam = bed.get_fam()
    
    if loadPhenos :
        # process phenotypes
        status= fam['status'] 
        
        if caseControl : y = recodeCaseControlQuantitative(status)  # as BCE now expects/outputs [N,1] shape, we must NOT use one-hot
        #if caseControl : # if trait is binary
        #    if recodeCaseControl : y = recodeCaseControlQuantitative(status) # if we want to treat it as a quantitative trait
        #    else : y = recodeCaseControlOneHot(status)# otherwise use the 'softmax' format
        else : 
            y = np.array(status)# otherwise just leave it alone
            y = y.reshape(-1,1) # enforce 2D , otherwise z-score wont work
    else : y = None # if there are no phenotypes
    
    # load genotypes from PLINK into RAM
    M = np.zeros( (bed.get_nb_samples(),bed.get_nb_markers()) , dtype =np.int8) # pre allocate a matrix with the correct size
    #loci_names = list()
    loci_names = np.empty(bed.get_nb_markers( ), dtype =   np.dtype((str, 25)) ) ## assume no rsid will be longer than 25 chars
    # Iterating over all loci
    counter = 0
    for loci_name, genotypes in bed:
        if replaceMissing : genotypes[genotypes==-1] = 0  # replace no calls with 0
        #if(M is None) : 
        #    M = genotypes
        #else : 
        #    M = np.column_stack( (M, genotypes) )
        M[:,counter] = genotypes # much faster to just paste this in (column stack is 10x slower)
        #loci_names.append(loci_name)
        loci_names[counter] = loci_name
        counter = counter +1
       
    # produce GCTA compatible ID list
    id_list = list()
    id_list.append( list(fam["fid"]))
    id_list.append( list(fam["iid"]))  
        
    return ( {"y" : y, "M" : M, "rsid" : loci_names.tolist(), "IDs" : id_list, "A1" : list(bim.a1), "A2" : list(bim.a2), "chrom" : list(bim.chrom), "pos" : list(bim.pos) } ) # return results

#A1_alleles = list(bim.a1)
def recodeCaseControlOneHot(status) : 
    # transcode phenotype into 1 hot array
    counter = 0
    y = np.zeros( (len(status),2)  )
    for i in range(len(status)) :
        result = np.array([0,0])
        result[ int(status[i] -1) ] =1 # transocde phenotype into one-of-k format: in plink case =2, control =1, so we set the element to 1 to at the inde-1 position
        y[i] =result  
        if(status[i] == 2) : counter = counter+1
    y = y.reshape(-1,2) # enforce 2D , otherwise z-score wont work
    return(y)


def recodeCaseControlQuantitative(status) : # transcode binary PLINK phenotype into a quantitative phenotype ( IE just turn the 1/2 into 0.0 and 1.0 )
    y = np.zeros( (len(status))  )
    for i in range(len(status)) :
        y[i] =(status[i] -1)  

    y = y.reshape(-1,1) # enforce 2D , otherwise z-score wont work
    return(y)



def recodeOneHotCaseControl(y) :
    # transcode 1 hot array into binary
    status = np.zeros( ( (len(y), 1) ) , dtype = int )
    for i in range(y.shape[0]) :
        
        status[ i ] =  int( np.argmax(y[i,:]) + 1 )
     
    
    return(status)


#location= "C:/!0datasets/adni/wgs/glmnettest/igap05_filtered_f1_train.pheno"
#results = loadPLINKPheno(location, caseControl = False, recodeCaseControl = False)
def loadPLINKPheno(location, caseControl = True, recodeCaseControl = True) :
    status = list()
    indisWithPhenos = list()
    with open(location, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            indisWithPhenos.append(itmp[1] )
            status.append( np.float32(itmp[2]) )
 
    if caseControl : y = recodeCaseControlQuantitative(status)  # as BCE now expects/outputs [N,1] shape, we must NOT use one-hot
        #if recodeCaseControl : y = recodeCaseControlQuantitative(status) # if we want to treat it as a quantitative trait
       # else : y = recodeCaseControlOneHot(status)# otherwise use the 'softmax' format                                      
    else : 
        y = np.array(status)# otherwise just leave it alone
        y = y.reshape(-1,1) # enforce 2D , otherwise z-score wont work
    
    return (  y.astype('float32') , indisWithPhenos )


# location = args.prs
# loads a PRS file from location, only those are kept which are found in plink file the lookups dictionary (based on chr, pos and A1/A2 matching)
def loadPRSFile(location, lookups, disableWeights = False) :   
    # go through text file which has signature ( with header):
    #   0         1          2                3              4
    # hm_chr	hm_pos	effect_allele	other_allele	effect_weight
    lookup_existing_dict = {}
    lookups_existing = list()
    betas = list()
    counter = 0
    numFlipped = 0
    with open(location, "r") as id:
        for i in id:
            counter+=1
            if(counter == 1) : continue # skip the header
            itmp = i.rstrip().split()

            chr = np.int32(itmp[0])
            pos = np.int32(itmp[1])
            lookup = tuple(np.array([chr, pos], np.int32))
            A1 = itmp[2]
            A2 = itmp[3]
            beta = np.float32(itmp[4])
            
            # check if it is present in the dictionary
            if lookup in lookups:
                 bims = lookups[lookup]
                 if bims[0] == A1 and bims[1] == A2 or bims[1] == A1 and bims[0] == A2: # check if it has the right alleles
                     if(A1 != bims[0]) : 
                         beta = -1 * beta #  flip allele effect if A1 wasn't the same between PRS and bim file
                         numFlipped+=1
                     lookups_existing.append(lookup) 
                     lookup_existing_dict[lookup] = len(lookups_existing) -1 # save the index of the SNP in the list
                     if disableWeights : betas.append(np.float32(1)) # if weights are disabled we just use 1
                     else : betas.append(beta)
    print("flipped number of Allele effects", numFlipped, "/", len(betas), flush=True)
    return  lookup_existing_dict, lookups_existing,  np.array(betas)
        

def loadPRSIndiFile(location, lookups) :  
    lookup_dict = {}
    for i in range(len(lookups) ):   
        lookup_dict[lookups[i] ] = i


    # go through text file which has signature ( with header):
    #   0         1          2          
    # IID         PHENO1     SCORE1_SUM
    scores = list()
    IDs = list()
    IDs_dict = {}
    counter = 0

    with open(location, "r") as id:
        for i in id:
            counter+=1
            if(counter == 1) : continue # skip the header
            itmp = i.rstrip().split()

            ID = itmp[0]
            score = np.float32(itmp[2])

            
            # check if it is present in the dictionary
            if ID in lookup_dict:
                IDs.append(ID)
                scores.append(score)
                IDs_dict[ID] = ID

    return  IDs, IDs_dict, np.array(scores)
    
 



# re
def validTrainSplit(validIDs, allIndis) :
    # first cache the valid indis list into a dict for fast lookups
    validIDs_dict = {}
    for i in range(len(validIDs) ):   
        validIDs_dict[validIDs[i] ] = i

    trainingSet = list()
    validSet = list()
    for i in range(len(allIndis) ): # loop the all indis array, to get indices referring to that
        lookup = allIndis[i]
        if lookup in validIDs_dict : validSet.append(i) # if indi in allindis is in the validation set, we save it
        else : trainingSet.append(i) # otherwise save it to the training set

    return(np.array(trainingSet, np.int32), np.array(validSet, np.int32))



            
# returns the indices in toBeMatchedList, that are found in the lookup_dict. Optionally, if indicesInOther, then it returns the indices in the other list (lookup_dict is assumed to have indices as values)
def matchIndices(toBeMatchedList, lookup_dict, indicesInOther = False) :
    matched_indices = list()
    for i in range(len(toBeMatchedList) ): # loop the array and not the dict for guaranteed order
        lookup = toBeMatchedList[i]
        if lookup in lookup_dict :
            if(indicesInOther ) : matched_indices.append(lookup_dict[lookup]) # if we want to get the indices in the other list
            else : matched_indices.append(i)
    return(np.array(matched_indices, np.int32))




# toBeMatchedList = IDs
# same as matchIndices, but compares against two dicts (but can only give indices according to the original input)
def matchIndices2(toBeMatchedList, lookup_dict, lookup_dict2) :
    matched_indices = list()
    #bothCriteriaMet = False
    for i in range(len(toBeMatchedList) ): # loop the array for guaranteed order
        lookup = toBeMatchedList[i]
        if lookup in lookup_dict and lookup in lookup_dict2:
            matched_indices.append(i)
    return(np.array(matched_indices, np.int32))


    # reorganise effects so that they are index-matched to the bim



# location=    args.covars_cont
 # location=    args.covars_factor
 # location=     args.covars_IDs
def loadNumericFromTextFile(location, forceIntDtype = False) :

    # go through whole thing to get the correct shape
    numLines = 0
    numCols = 0;
    with open(location, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            numLines+=1
            numCols = len(itmp)

    dataType = np.float32
    if forceIntDtype : dataType = np.int32
    data = np.zeros( (numLines,numCols) , dtype =dataType) # pre allocate a matrix with the correct size
    #data.shape
   
    with open(location, "r") as id:
        counter = 0;
        for i in id:
            itmp = i.rstrip().split()
            status = list()

            for j in range(len(itmp) ):
                status.append( np.float32(itmp[j]) )

            data[counter,:] =  np.array(status)
            counter+=1

    return ( data )

# location = args.covars_IDs
def loadTextList(location) :
    status = list()
    with open(location, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            status.append(itmp[0] )
    return (  status )


def loadIntegerList(location) :
    status = list()
    with open(location, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            status.append( np.int32( itmp[0]) )
    return (  status )
    


# column binds Genotype matrices (IE stacks them next to each other)
#M_list = [ results2["M"], results["M"] ]
#allM = concatChroms(M_list)
def concatChroms(M_list) : 
    return(np.concatenate(M_list, axis=1))
    
##################################################################################################
# Generic data IO (IE any kind of matrix or array)
##################################################################################################

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
    

# loads matrix from disk ( that was written by the above)
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


# writes an array onto disk ( 2 files, 1 text file that stores the length, the other the binary)
def writeVectorToDisk(location,data, dataType ="float32") :
    d_code = get_d_code(dataType)
    
    # write the dimensions onto disk
    with open(location + ".id", "w") as idFile: 
        idFile.write( str( len(data) ) )

    flatData = struct.pack(d_code*len(data),*data  )
    with open(location + ".bin", "wb") as flat_File: 
        flat_File.write(flatData) 
    

    
# loads array from disk ( that was written by the above function)
def loadVectorFromDisk(location, dataType ="float32") :
    d_code = get_d_code(dataType)
    
    # load id file to get dimensions
    with open(location + ".id", "r") as idFile:
        itmp = idFile.readline().rstrip()
        totalNum = int(itmp) # how many elements to expect in the binary in total

    # open binary file
    with open(location + ".bin", "rb") as BinFile:
        BinFileContent = BinFile.read()
    
    # reformat data into correct dimensions
    flat = np.array( struct.unpack(d_code*totalNum, BinFileContent  ), dtype = dataType  )
    return(flat)
    
    
def loadIndices(location) :
    fileName = location
    indices = list()
    with open(fileName, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            indices.append( np.int(itmp[0]) )
      
    return ( np.array(indices) )
    
def loadsnpIDs(location) :
    fileName = location
    snpIDs = list()
    with open(fileName, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            snpIDs.append( itmp[0] )
      
    return ( np.array(snpIDs) )

def writeSNPeffects(location,SNPIds, SNP_effects) :
    fileName = location
    with open(fileName, "w") as file: 
        for i in range( len(SNPIds) ):
            file.write( str(SNPIds[i]) + "\t" + str(SNP_effects[i])  + "\n")   

##################################################################################################
# GCTA formatted kinship, as these are Lower Triangle matrices they have to be treated separately
##################################################################################################


# writes a kinship matrix onto disk in the GCTA format
def writeGCTA_GRM(location,K, id_list, numSNPs) :
    # generate filenames
    BinFileName = location+".grm.bin"
    IDFileName = location+".grm.id"
    NFileName = location+".grm.N.bin"
    
    # how many to write
    numIndis = len(id_list[0])
    numPairs = int(  numIndis*(numIndis+1)/2 )
    N_vals = [numSNPs] * numPairs # the number of surviving SNPs used, the number of Individuals 
    
    # write N file
    NFileData = struct.pack("f"*numPairs,*N_vals  )
    del N_vals; gc.collect();
    with open(NFileName, "wb") as NFile: 
        NFile.write(NFileData) 
    del NFileData; gc.collect();
             
                                         
    # write ID list
    with open(IDFileName, "w") as idFile: 
        for i in range( len(id_list[0]) ):
            idFile.write(id_list[0][i] + "\t" + id_list[1][i] + "\n") 
    gc.collect();

    
    # write GRM to disk: unravel the GRM's lower triangle
    grm_indices = np.tril_indices(numIndis)
    grm_pairs= list (K[grm_indices[0], grm_indices[1]] )
    del grm_indices; gc.collect();
    
    # write out binary
    GRMFileData = struct.pack("f"*numPairs,*grm_pairs  )
    del grm_pairs; gc.collect();
    with open(BinFileName, "wb") as GRM_File: 
        GRM_File.write(GRMFileData) 
    

# loads a kinship matrix from disk in the GCTA format ( written out by the above function)
def loadGCTA_GRM(location) :
    # generate filenames
    BinFileName = location+".grm.bin"
    IDFileName = location+".grm.id"
    NFileName = location+".grm.N.bin"
    
    # load the familiy/individual IDs
    id_list = list()
    id_list.append([])
    id_list.append([])
    with open(IDFileName, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            id_list[0].append(itmp[0])
            id_list[1].append(itmp[1])
    id.close()
    
    # how many peopple in total, IE the number of entries in the binary file
    n_subj = len(id_list[0])
    nn = int(  n_subj*(n_subj+1)/2 )
    
    # load Binary
    with open(BinFileName, "rb") as BinFile:
        BinFileContent = BinFile.read()
    BinFile.close()
    
    # reformat it into the correct dimensions ( as only the lower triangle was stored flattened)
    K = np.zeros((n_subj, n_subj))
    
    grm_vals = list(struct.unpack("f"*nn, BinFileContent  ))
    inds = np.tril_indices_from(K)
    K[inds] = grm_vals # reconstruct the full matrix from the LT
    K[(inds[1], inds[0])] = grm_vals
      
    # load the 'N' file
    with open(NFileName, "rb") as NFile:
        NFileContent = NFile.read()
    NFile.close()
    
    N_vals = list(struct.unpack("f"*nn, NFileContent  ))
    
    gc.collect();
    return ( {"K" : K, "ids" : id_list, "N" :N_vals} ) # return results


def load_N_fromGCTA_GRM(location) :

    NFileName = location+".grm.N.bin"

    # load the 'N' file
    with open(NFileName, "rb") as NFile:
        NFileContent = NFile.read(4) # only load the first 4 bytes = 16 bit float, as unlike GCTA I use the same number of SNPs for all
        
    N = struct.unpack("f", NFileContent  )[0]
   
    return ( N ) 


# writes human readable text file of a Variance Compnent Analysis results
def writeVCResults(location,allModels ) :

    # write ID list
    with open(location, "w") as targetFile: 
        
        for i in range( len(allModels) ):
            line = "Model with "+ str( len(allModels[i]["vc"]) ) + " genetic variance component(s):"
            print(line) # want to save output to console too..
            targetFile.write(line + "\n")
            
            line = "BIC: " + str( allModels[i]["bic"] )
            print(line)
            targetFile.write(line + "\n") 
            
            line = "h2: " + str( allModels[i]["h2"] )
            print(line)
            targetFile.write(line + "\n") 
            
            line = "h2_CI95 lb: " +  str( allModels[i]["h2_lb"] ) + " / ub: " +  str( allModels[i]["h2_ub"] )
            print(line)
            targetFile.write(line + "\n")
            
            line = "Ve: " +  str( allModels[i]["ve"] )
            print(line)
            targetFile.write(line + "\n")
            
            epistasisCounter = 2
            for j in range( len(allModels[i]["vc"]) ) : # go through all variance components
                note = ""
                if j == 0 : note ="(additive)" # first element is always the additive VC
                elif j == allModels[i]["domcol"] : note ="(dominance)" # if j is the same col as the dominance effect
                else : 
                    note ="(epistasis "+str(epistasisCounter)+"-way)"
                    epistasisCounter = epistasisCounter +1
                    
                line = "VC"+str(j+1)+" "+note+": "+ str( allModels[i]["vc"][j] ) + " / p-value: " +  str( allModels[i]["p"][j] ) + " / vc_se: " +  str( allModels[i]["vc_se"][j] )
                print(line + "\n")
                targetFile.write(line + "\n")

            line = "_________________________________" 
            print(line)
            targetFile.write(line + "\n") 

##################################################################################################
# write regions
##################################################################################################



def writeRegionData(location,regions, deltas) :
    fileName = location
    with open(fileName, "w") as file: 
        for i in range( len(deltas) ):
            file.write( str(regions[i][0]) + "\t" + str(regions[i][1]) + "\t" + str(deltas[i])  + "\n")   

def loadRegionData(location) :
    fileName = location
    regions = list()
    deltas = list()
    
    with open(fileName, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            regions.append( np.array( [int(itmp[0]) , int(itmp[1])] ) )
            deltas.append( np.float64(itmp[2]) )
            
    return ( {"REGIONS":regions, "DELTAS":deltas } )

##################################################################################################
# Summary Stats
##################################################################################################

#location ="C:/!0datasets/adni/wgs/stage1_p05_short_ordered.txt"
def loadSummaryStats(location) :
    fileName = location
    summaryData = list()

    with open(fileName, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            summaryData.append(  np.float64(itmp[2]) ) # the Betas are stored in the 3rd column
 
    # process them into priors
    # take their abs ( as it doesnt matter if its negative or positive)
    summaryData = np.abs(summaryData)
    
    # scale them to be 0-1
    max_x = np.max(summaryData)
    min_x = np.min(summaryData)
    denominator = (max_x  - min_x)
    if denominator != 0 : summaryData = (summaryData -min_x )/ denominator # watch out fo division by 0
   
    return ( summaryData )

#import numpy as np
#location ="C:/!0datasets/adni/wgs/h2/stage1score05_filtered_001_ordered_pruned.txt"
#priors = summaryData
def loadSummaryStats_noscale(location) :
    fileName = location
    summaryData = list()

    with open(fileName, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            summaryData.append(  np.float64(itmp[2]) ) # the Betas are stored in the 3rd column

    return ( summaryData )


##################################################################################################
# LDAK helper
##################################################################################################

def loadLDAKWeights(location) :
    fileName = location
    weights = list()
    
    counter = 0
    with open(fileName, "r") as id:
        for i in id:
            if counter > 0 : # weight file has header
                itmp = i.rstrip().split()
                weights.append( np.float64(itmp[1]) )
                
            counter = counter +1
            
    return ( weights )

def loadKinshipLocations(location) :
    fileName = location
    kinshipLocations = list()

    with open(fileName, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            kinshipLocations.append(itmp[0] )
 
            
    return ( kinshipLocations )


# loads the file that stores how many regions we have, and then loads each file in turn and adds their rsids into
def loadLDAKRegions(location) :
    fileName = location + "region_number"
    
    numRegions = -1
    with open(fileName, "r") as NFile:
        #itmp = NFile.rstrip().split()
        numRegions = np.int(NFile.readline())

    numRegions = numRegions +1 #  as we have the background region, which is 'region0'
    
    allRegions = list() # a list of lists
    for j in range(numRegions) :
        fileName = location + "region" + str(j)
     
        nextRegion = list()
        allRegions.append(nextRegion)
        with open(fileName, "r") as region:
            for i in region:
                nextRegion.append(i.rstrip().split()[0]) # each line in the file is the rsid

                
    return ( allRegions )


def loadLDAKRegionsDeltas(location) :# takes an reml file
    MAXDELTA = np.exp(10)
    allDeltas = list()
    with open(location, "r") as remlResults:
        for i in remlResults:
            if i.find("Her_") > -1 : # 
                itmp = i.rstrip().split() 
                h2 = np.float64( itmp[1]) # the 2nd item is the h2
                if h2 == 0 : delta  = MAXDELTA
                else : delta =   (1-h2) /h2
                allDeltas.append(delta)
                # h2 = Va/Vphe  ~ = Va/1   = Va 
                # Vphe = Va + Ve -> Ve = Vphe - Va  = Ve = 1 -h2
                # delta = Ve/Va == (1-Va) / Va  == (1-h2) / h2
                
    return(allDeltas)


# lines that contain the h2 results all contain "Her_"
# Her_K1 0.000000 0.000000
# Her_R1 0.000135 0.000773




    
##################################################################################################
# write/load eigen Summaries of kinships
##################################################################################################

    
    
def writeEigenSum(location, eigSums) :
    for i in range( len(eigSums) ):
        writeMatrixToDisk(location + "eigVec" + str(i+1) , eigSums[i].vectors)
        writeVectorToDisk(location + "eigVal" + str(i+1), eigSums[i].values)


 
def loadEigenSum(location) :
    eigSums = list()
    moreFiles = True
    counter = 1
    while moreFiles : # keep going until we run out of files to load
        currentLocation = location + "eigVec" + str(counter) + ".id"  # files are expected to be named as eigVec1.id/bin etc
        my_file = Path(currentLocation)
        if my_file.is_file(): # check if it exists
        
            filename = location + "eigVec" + str(counter)
            vectors = loadMatrixFromDisk(filename)
            filename = location + "eigVal" + str(counter)
            values = loadVectorFromDisk(filename)
            
            eigenSummary = collections.namedtuple('values', 'vectors')
            eigenSummary.values = values
            eigenSummary.vectors =vectors
 
            eigSums.append( eigenSummary ) 
        else : moreFiles = False
        counter = counter +1
        
    return(eigSums)
    
    

    
    
    
    
##################################################################################################
##################################################################################################

# location="dummyregions" # local test
#deltas = [1,2]
#regions = list()
#regions.append( np.array([0,50]))
#regions.append( np.array([25,75]))
#writeRegionData(location,regions, deltas)
# results = loadRegionData(location)

# local testing
#gcta_data2 = loadGCTA_GRM(location)
#location = "data/gcta/output/wtccc2_hg19_toy"
# writeGCTA_GRM(location, gcta_data["K"], gcta_data["ids"], gcta_data["N"])
    

    
# local test
#location = "data/testMatrix"
#data = np.random.normal(size=[ 100,100], scale=1)

#location = "data/testVec"
#dataVec = np.random.normal( size=50, scale=1)