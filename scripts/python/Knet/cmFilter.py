# -*- coding: utf-8 -*-
"""
Created on Sun May  5 18:43:39 2019

@author: mk23
"""

# Given a UKBB bgen variant stats file, it gives back a range that equals a 500SNP flanking window around 2 variants

import argparse
import numpy as np

from itertools import combinations



from pathlib import Path



def findCMBlock(BpA,boundaries) : 
    for i in range(len(boundaries)) :
        if (BpA >= boundaries[i][0] and  BpA <= boundaries[i][1] ) :
            return(i)
    
    return(-1)

##################################################################

def processPRS(cmfile, sumstats,out,header,bpThreshold) :

    # load centi Morgan LD recombination blocks file
    chroms=list()
    lastChrom = -1
    currentChrom = -1
    counter = 0
    with open(cmfile, "r") as id:
        for i in id:
          #  print (i)
            counter+=1
            if (counter != 1) : # there is header
                itmp = i.rstrip().split()
                chrom=int(itmp[0])
                cm_start=int(itmp[1])
                cm_end=int(itmp[2])
                
                if(chrom != lastChrom) :
                    #print("found new chrom:", chrom)
                    currentChrom = []
                    chroms.append(currentChrom)
                    lastChrom = chrom
                    
                    
                currentChrom.append([cm_start,cm_end])
    
    counter = 0
    sameCMBPredictors = []
    # parse the summary stats file, it is assumed 
    with open(sumstats, "r") as id:
        for i in id:
            #print (i)

            counter+=1
            print("checking ", counter)
            if (counter != 1 or header == False) : # there is header
                itmp = i.rstrip().split()
                chromA=int(itmp[2])  
                chromB=int(itmp[3])   
                
                if(chromA == chromB) : # if they are on the same chrom (this should always be the case if we got a pre-filtered assoc file)
                    #print("same chrom")
                    # get their bp positions

                    BpA=int(itmp[4])  
                    BpB=int(itmp[5])  
                    
                    # check if they are within the same recombination block
                    
                    # find the cMs for this chrom
                    boundaries = chroms[(chromA-1)] # -1 as it is 0 indexed
                    predictor1_block = findCMBlock(BpA,boundaries)
                    predictor2_block = findCMBlock(BpB,boundaries)
                    
                    distance = abs(BpA-BpB)
                    
                    # bpThreshold if we they are closer than this then they are automatically considered too close ( it is possible that )                
                    # examples when this happened
                    # chrom 11abs(67372946-67636004) = 263058
                    # PWAS: 8052	162, chrom 19: abs(55993507-56015159) = 21652  <- they are like 20K Bases away but right at the boundaries, so did not get picked up by the system

                 
                    
                    if (predictor1_block == predictor2_block or distance <= bpThreshold) :
                        predictor1=itmp[6]
                        predictor2=itmp[7]
                        #print(predictor1, " and ", predictor2 , "are in the same cM block, we exclude them")
                        sameCMBPredictors.append([predictor1,predictor2])
                        
                        
                        
    with open(out, "w") as resultsFile: 
        #resultsFile.write("predictor1" + "\t" + "predictor2"  + "\n" )
        for i in range(len(sameCMBPredictors)  ) :
            resultsFile.write(str(sameCMBPredictors[i][0]) + "\t" + str(sameCMBPredictors[i][1]) + "\n" )

    print("written ", len(sameCMBPredictors) ," pairs of predictors that are in the same cM block to: ", out)

#myString="CACAGACGCAGGGCCCACCACTACCTTGAGGACGGGGCCAAGGAAAGTGTCCTTCTTGTGCCGGCAGACTCTGCTGAGGACCAGGAAAGCCAGCACCCGCAGAGACTCTTCCCCAGTGCTCCATACGAT"
#myString="AAGAGCAAGCACATTACGAGCCACACTTGGTGCTGGTGCTGTGGCATGCAACACTACCTAATGCGAGAGAAAGATGTGAGCAAATATCCGGAATATACAAATATAATACAAAAATACAAATATAATCA"
#49431 / 101476
#len(myString)
#-6.977e+02
# 154+ 647
#start=1000-64
#end=1000+64
#end-start
#1028- 900
# header=False     
# cmfile='C:/softwares/Cluster/GIANT/miniPRS/epistasis/all.1cM.tab'
#out='C:/softwares/Cluster/GIANT/miniPRS/epistasis/epistasis_results_ALL_samechrom_short_fail'         
#sumstats='C:/softwares/Cluster/GIANT/miniPRS/epistasis/epistasis_results_ALL_samechrom_short'
def runMain(args) :
    print('Processing sumstats from ',  args.sumstats)    
    out = args.out
    cmfile = args.cmfile
    sumstats = args.sumstats
    header= False
    if  args.header == "1" : header = True
    bpThreshold = args.bpThreshold * 1000
    processPRS(cmfile, sumstats,out,header,bpThreshold)
  #  (675 - 628) /2

if __name__ == '__main__':   
    parser = argparse.ArgumentParser()

    parser.add_argument("--out",required=True, help='an output location is always required')
    parser.add_argument("--cmfile",required=True, help='Path to the cM map of the genome, signature: chr   start   end')  
    parser.add_argument("--sumstats",required=True, help='Path to the epistasis association file, produced by OLS_epistasis (NO HEADER)')  
    parser.add_argument("--header",required=False, help='If the summary stats file is assumed to have a heder (default 0: No header)', default="0")  
    parser.add_argument("--bpThreshold",required=False, type=int, help='distance within markers are automatically considered to be too close (default 500 kb)', default=500)  

    parser.set_defaults(func=runMain)
        
    args = parser.parse_args()
    args.func(args)








