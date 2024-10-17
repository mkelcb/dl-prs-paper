# -*- coding: utf-8 -*-
"""

given a list of pairs of variants whose combination failed QC
- it removes the least number of variants so that all combinations go away
@author: mk23
"""

# Given a UKBB bgen variant stats file, it gives back a range that equals a 500SNP flanking window around 2 variants

import argparse
import numpy as np

from itertools import combinations



from pathlib import Path



# loads pairs which failed QC, signature will be like
# rs1234    rs5678
def loadFailedVars(failedPairsLoc) :
    faledPairs=list()

  #  counter = 0
    with open(failedPairsLoc, "r") as id:
        for i in id:
            itmp = i.rstrip().split()
            faledPairs.append( [itmp[0],itmp[1]] )

          #  counter += 1
 
    return(faledPairs)


##################################################################
# A    B
# B    C
# B    D
# B    E
# - it is enough to exclude just B, to remove all bad combinations
# - algo to remove only the least number of predictors to remove all bad combinations:
# - start from top, remove one arbitrarily, then loop, rest and remove all combinations that included that 
# - then move onto the next one

#
#failedPairs = [ ["A","B"],["B","C"],["B","D"],["B","E"] ]

#failedPairs = [ ["A","B"],["B","C"],["B","D"],["B","E"],["A","F"] ]

#failedPairs = [ ["A","B"],["B","C"],["B","D"],["B","D"],["B","E"],["A","F"] ]
#failedPairs = [ ["A","B"],["B","C"],["B","D"],["B","F"],["B","E"],["A","F"] ]
def process(failedPairs,out) :

    predictorsToRemove = [] # the final list of predictors that if removed all bad combinations will go away
    while(len(failedPairs) > 0) : # keep going until we run out
        currentPredictorPair = failedPairs.pop(0) # remove the first

        currentPredictor1 = currentPredictorPair[0]
        currentPredictor2 = currentPredictorPair[1]

        # we don't know which of the pair is featured in more pairs, 
        # so we check both to remove the one that is featured in more pairs
        print("going to remove " , str(currentPredictor1), "and/or", str(currentPredictor2), " num pairs to check left ", len(failedPairs ) )
        
        nextPredictorsToRemove_ifRemoveFirst = []
        for i in range(len(failedPairs)) : # go through the remaining pairs, 
            nextPair = failedPairs[i]
            if nextPair[0] == currentPredictor1  or nextPair[1] == currentPredictor1 : # and if we find any pair where the currentPredictor to be removed is a member, we can then remove that pair from the list to be checked
                nextPredictorsToRemove_ifRemoveFirst.append(i) # add the predictors to be removed to the list
                
        nextPredictorsToRemove_ifRemoveSecond = []
        for i in range(len(failedPairs)) : # go through the remaining pairs, 
            nextPair = failedPairs[i]
            if nextPair[0] == currentPredictor2  or nextPair[1] == currentPredictor2 : # and if we find any pair where the currentPredictor to be removed is a member, we can then remove that pair from the list to be checked
                if i not in nextPredictorsToRemove_ifRemoveFirst: # dont try to remove the same element twice (should not happen)
                    nextPredictorsToRemove_ifRemoveSecond.append(i) # add the predictors to be removed to the list
               
        
        # potentially remove one or both predictors, in case they were featured in an additional pair
        if len(nextPredictorsToRemove_ifRemoveFirst) > 0 :
            predictorsToRemove.append(currentPredictor1)  
            print("removing",currentPredictor1 ," as it was in num pairs: ", len(nextPredictorsToRemove_ifRemoveFirst) )

        
        if len(nextPredictorsToRemove_ifRemoveSecond) > 0 :
            predictorsToRemove.append(currentPredictor2)     
            print("removing",currentPredictor2 ," as it was in num pairs: ", len(nextPredictorsToRemove_ifRemoveSecond) )
        
       # nextPredictorsToRemove_ifRemoveFirst = [1,2,5,0]
       # nextPredictorsToRemove_ifRemoveSecond = [3,6,4]
        
       # join the two lists, then sort them descending so that we start with the last index
        nextPredictorsToRemove_ifRemoveFirst.extend(nextPredictorsToRemove_ifRemoveSecond)
        nextPredictorsToRemove_ifRemoveFirst.sort()
        nextPredictorsToRemove_ifRemoveFirst.reverse()
          
        # remove them from the 'list left to do' by index, this is safe, as we always remove the last index first
        while(len(nextPredictorsToRemove_ifRemoveFirst) > 0 and len(failedPairs) > 0) : # go through list and remove all pairs
            PairToremove=nextPredictorsToRemove_ifRemoveFirst.pop(0)
            print("removing index",PairToremove )
            failedPairs.pop(PairToremove)
 
            
            
    with open(out, "w") as resultsFile: 
        #resultsFile.write("predictor1" + "\t" + "predictor2"  + "\n" )
        for i in range(len(predictorsToRemove)  ) :
            resultsFile.write(str(predictorsToRemove[i])+ "\n" )

    print("written ", len(predictorsToRemove) ," predictors needed to be removed to: ", out)


# header=False     
# cmfile='C:/softwares/Cluster/GIANT/miniPRS/epistasis/all.1cM.tab'
#out='C:/softwares/Cluster/GIANT/miniPRS/predsToRemove'         
#sumstats='C:/softwares/Cluster/GIANT/miniPRS/epistasis/epistasis_results_ALL_samechrom_short'
def runMain(args) :
    print('Processing failedPairsLoc from ',  args.failedPairsLoc)    
    out = args.out
    failedPairsLoc = args.failedPairsLoc

    failedPairs = loadFailedVars(failedPairsLoc)
    process(failedPairs,out)


if __name__ == '__main__':   
    parser = argparse.ArgumentParser()

    parser.add_argument("--out",required=True, help='an output location is always required')
    parser.add_argument("--failedPairsLoc",required=True, help='Path to the QC-failed combinations of predictors, signature: rs1234    rs5678')  

    parser.set_defaults(func=runMain)
        
    args = parser.parse_args()
    args.func(args)



