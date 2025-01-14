# fits a regression based linear models for the PRS GXE project: Cov-only and PRS+Cov
#options(error=traceback)
# grab command line arguments
# options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 9) {stop("not enough Arguments received")} 


trainValidPhenosLoc= args[1] # the location for all the phenos, that has data for the trainvalid indis
testIndisLoc = args[2] # the location for the file for the test indis
validIndisLoc = args[3] # the location for the file for the valid indis
PRSloc= args[4] # location for the file for all of the PRS
CovIDsLoc= args[5] # location for the file for the IDs for the covariates
CovNumericsLoc= args[6] # location for the numerical covariates
CovFactorsLoc= args[7] # location for the factor covariates
isBin= args[8] == "1" # flag to decide if binary (1) or quantitative pheno (0)
outputLoc = args[9] # where to write results

# # DEBUG VARS
# trainValidPhenosLoc= "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/MONDO_0002009_all_Sex_trainValid"
# testIndisLoc = "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/test_keep_largecases2"
# validIndisLoc = "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/valid_keep_largecases2"
# PRSloc= "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/MONDO_0002009/PRS_all.sscore"
# CovIDsLoc= "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_ids"
# CovNumericsLoc= "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_numericBinary_all"
# CovFactorsLoc= "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/data/raw/COVS_factors_all"
# isBin= T
# outputLoc = "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/results/MONDO_0002009res_reg"


# Load data:
trainValidPhenos=  read.table(trainValidPhenosLoc ,header=F) 
testIndis = read.table(testIndisLoc ,header=F) 
validIndis = read.table(validIndisLoc ,header=F) 
PRS= read.table(PRSloc ,header=T)  # IID     PHENO1  SCORE1_SUM
CovIDs=   read.table(CovIDsLoc ,header=F)
CovNumerics= read.table(CovNumericsLoc ,header=F) 
CovFactors= read.table(CovFactorsLoc ,header=F)
# sapply(CovFactors, class) # this looks numeric, but must be cast to factor
for(i in 1:ncol(CovFactors)) {
  CovFactors[,i] <- as.factor(CovFactors[,i]) 
}
factorColnames = paste0("fact", 1:ncol(CovFactors), "_level") # add the '_level' postfix, so that when we want to find significant predictors, we have a way of splitting between the stem and the factor level generated by R
colnames(CovFactors) = factorColnames # need unique names for these, so that we could later use those to exclude them

# add column for IDs for Cov numeric/factors
colnames(CovIDs) = c("IID")
CovNumerics = cbind.data.frame(CovIDs,CovNumerics)
CovFactors = cbind.data.frame(CovIDs,CovFactors)

# remove from trainValid the valid
# find the indices of the test indis
ValidIndices = match(validIndis$V1, trainValidPhenos$V1) # find the indices in the trainvalid data that refer to the valid set
ValidIndices= ValidIndices[is.na(ValidIndices) == F] # drop NAs (indis that could not be found)
#trainValidIndices = 1:nrow(trainValidPhenos)

# remove from trainValid the valid
trainIndis = trainValidPhenos[-ValidIndices,]

# create design matrix: combine PRS+Cov-numeric and Cov-factor (must keep this 'as factor')
# first create ALL of the phenotypes, ie the 125K train, test,
allIndis = c(trainValidPhenos$V1,testIndis$V1)
allIndis = as.data.frame(allIndis)

# get the pheno and the PRS : $PHENO1 ,  $SCORE1_SUM
PRS_pheno = merge(PRS, allIndis, by.x = "IID" , by.y ="allIndis" )

PRS_pheno_covs = merge(PRS_pheno, CovNumerics, by.x = "IID" , by.y ="IID" )
PRS_pheno_covs = merge(PRS_pheno_covs, CovFactors, by.x = "IID" , by.y ="IID" )
nrow(PRS_pheno_covs) # 124996, we lost a few people??
ncol(PRS_pheno_covs) # 201
PRS_pheno_covs = as.data.frame(PRS_pheno_covs)


# subset to test set with supplied data (PRS/Covs)
testIndices = match(testIndis$V1, PRS_pheno_covs$IID) # find the indices in the targetDF that refer to the test set
testIndices= testIndices[is.na(testIndices) == F] # drop NAs (indis that could not be found)
test_pheno = PRS_pheno_covs$PHENO1[testIndices]
test_IDs = PRS_pheno_covs$IID[testIndices]
test_PRS = PRS_pheno_covs$SCORE1_SUM[testIndices]

# subset to train set with supplied data (PRS/Covs)
trainIndices = match(trainIndis$V1, PRS_pheno_covs$IID) # find the indices in the targetDF that refer to the test set
trainIndices= trainIndices[is.na(trainIndices) == F] # drop NAs (indis that could not be found)
train_pheno = PRS_pheno_covs$PHENO1[trainIndices]
#train_IDs = PRS_pheno_covs$IID[trainIndices]

## drop unused levels, otherwise will get "contrasts can be applied only to factors with 2 or more levels" error
# only do this AFTER we have subset it to the actual used indis, otherwise we may get some that only loose factor levels in the subset
CovFactors_used = PRS_pheno_covs[,factorColnames]

# subset this into the 3 different splits, and make sure that none of them violates the above:
#dataSplit=test_data_dummy
getValidFactorNames= function(dataSplit) {
  dim(dataSplit)
  dataSplit <- lapply(dataSplit, droplevels)
  dataSplit = as.data.frame(dataSplit)
  indices_to_remove = sapply(dataSplit, nlevels) <= 1    # remove factors with 1 level # https://stackoverflow.com/questions/17995195/remove-variables-with-factor-level-1
  indices_to_remove = which(indices_to_remove == 1)
  dataSplit = dataSplit[,-indices_to_remove]
  dim(dataSplit)
  return(colnames(dataSplit))
}

notValidIndices = c(trainIndices,testIndices)
test_factors_dummy = CovFactors_used[testIndices,]
train_factors_dummy = CovFactors_used[trainIndices,]
valid_factors_dummy = CovFactors_used[-notValidIndices,]

test_factors = getValidFactorNames(test_factors_dummy)
train_factors = getValidFactorNames(test_factors_dummy)
valid_factors = getValidFactorNames(valid_factors_dummy)
# intersect the 3 list of factors to find out which were good in all 3
test_train_factors = intersect(test_factors,train_factors) 
test_train_valid_factors = intersect(test_train_factors,valid_factors) 


badFactors = setdiff( factorColnames, test_train_valid_factors)
badFactors
#PRS_pheno_covs_orig = PRS_pheno_covs
# find out which factors were bad that we had to remove them
if(length(badFactors) > 0 ) {
  badFactors_indices = match(badFactors, colnames(PRS_pheno_covs))
  PRS_pheno_covs = PRS_pheno_covs[,-badFactors_indices]
  # now also remove them from the 3 subsets
  badFactors_indices = match(badFactors, colnames(test_factors_dummy))
  test_factors_dummy = test_factors_dummy[,-badFactors_indices]
  train_factors_dummy = train_factors_dummy[,-badFactors_indices]
  valid_factors_dummy = valid_factors_dummy[,-badFactors_indices]
  
}

# still need to cater for the situation where we kept a factor, as it has more than 1 level in each split, but these levels are NOT the same
# training set: "a" "b" "c"
# test set: "a", "c","d"
#i=107
colnam = colnames(train_factors_dummy)
for (i in 1:length(train_factors_dummy)) {
  
  train_table = table(train_factors_dummy[,colnam[i]])
  test_table = table(test_factors_dummy[,colnam[i]])
  valid_table = table(valid_factors_dummy[,colnam[i]])
  all_table = table(PRS_pheno_covs[,colnam[i]])
  
  # find any unused levels in either train/test
  unused_train = as.numeric(which(train_table == 0))
  unused_test = as.numeric(which(test_table == 0))
  unused_valid = as.numeric(which(valid_table == 0))
  unused_train_test = union(unused_train, unused_test)
  unused_train_test_valid = union(unused_train_test, unused_valid)
  # if there were any unused factors, we set these to be the mode in the 'all' and then drop levels
  if ( length(unused_train_test_valid) ) { 
    unusedFactors = unlist(dimnames(all_table))[unused_train_test_valid]
    modevalue = names(sort(-table(PRS_pheno_covs[,colnam[i]])))[1]
    print(paste0("removing from factor ", colnam[i], " unused factor levels ", unusedFactors, " replacing them by mode value: ", modevalue))
    levels(PRS_pheno_covs[,colnam[i]])[levels(PRS_pheno_covs[,colnam[i]]) == unusedFactors] <- modevalue 
    PRS_pheno_covs[,colnam[i]] <- droplevels(PRS_pheno_covs[,colnam[i]])
    
    # what happens if after the mode is set, we only have 1 level: we should just remove the whole column then
    if ( nlevels(PRS_pheno_covs[,colnam[i]]) <= 1) {
      indexToRemove = match(colnam[i], colnames(PRS_pheno_covs))
      PRS_pheno_covs = PRS_pheno_covs[, -indexToRemove]
      print(paste0("As factor ", colnam[i], " had only 1 level left, we remove it completely "))
    }
  }
}


# create the final design matrices:
test_data =PRS_pheno_covs[testIndices,3:ncol(PRS_pheno_covs)]
train_data =PRS_pheno_covs[trainIndices,3:ncol(PRS_pheno_covs)]



fitModel = function(train_data, test_data, scenar, CovarWriteOut = F) {
  
# fit regression model on train set
# Covs + PRS
if(isBin) {
  train_pheno_01 = train_pheno -1 # because PLINK has cases as 2, and controls as 1, but R does not like that
  model = glm(train_pheno_01 ~ .,family=binomial(link='logit'),data=train_data)
} else {
  model = lm(train_pheno ~ .,data=train_data)
}
#summary(model)

# fit prediction model on test set, write out in format: "ID\tpheno\tPRS" \n 1000101 1       0.120438"
y_hat = predict(model, newdata=test_data)
#print( paste0(scenar," cor: ", cor(y_hat, test_pheno) ) ) #  "ugaaa cor: 0.160021082132457"
origRes = cbind(y_hat, test_pheno)
origRes = na.omit(origRes)
print( paste0(scenar," cor: ", cor(origRes[,1], origRes[,2])) ) #  "ugaaa cor: 0.160021082132457"


#length( which(is.na(test_pheno) ) )
#length( which(is.na(y_hat) ) )
# y_hat_orig = y_hat

# write results to disk
testResults = cbind.data.frame(test_IDs, test_pheno, y_hat	)
colnames(testResults) = c("ID","pheno","PRS")
outL = paste0(outputLoc, scenar)
write.table(testResults,file=outL, quote = F, row.names = F, sep = "\t",append=F)
print(paste("written: results to", outL ))

# also write the final covariates data used, so that NNs only need to work on non-redundant data:
if(CovarWriteOut) {
  # extract the model's coefficients
  coefs <- summary(model)$coefficients
  coefs = coefs[2:nrow(coefs),] # 2: as we do NOT want the intercept no matter what
  # Identify significant variables
  sig_predictors <- rownames(coefs)[which(coefs[, 4] < 0.05)] 
  
  # because multi-level factors get a new dummy variable for each level (ie the 3 factors of "a", "b" and "c" will be turned into 3 different predictors)
  # to find out which of the original predictors were used, we need to remove the factor levels to recover original column names
  sig_predictors_cleaned = c()
  for(i in 1:length(sig_predictors)) {
    #print(paste0(sig_predictors[i], " is factor: ",grepl("fact", sig_predictors[i])  ) )
    if(grepl("fact", sig_predictors[i])) { # if its a factor, then it will have the "fact" in it, as I've renamed all column names like that
      
      # split the generated dummy variable name into its components, and get the first part, which will be the original column name
      colname = unlist (strsplit(sig_predictors[i],"_level"))[1]
      colname = paste0(colname, "_level") # add this back, otherwise we would not be able to find the original column names
      # now check if we don't already have this, and if not, then we add it
      if( length(sig_predictors_cleaned) == 0 || colname %in% sig_predictors_cleaned == FALSE) {
        sig_predictors_cleaned= c(sig_predictors_cleaned,colname )
        #print(paste0(sig_predictors[i], " is factor: ",grepl("fact", sig_predictors[i])  , " / colname:", colname) )
      }
      
    } else { # if its not a factor we keep it anyway
      sig_predictors_cleaned = c(sig_predictors_cleaned, sig_predictors[i])
    }
  }
  
  # write out the IDs separately
  write.table(PRS_pheno_covs$IID,file=paste0(outL,"_","IDs"), col.names = F, quote = F, row.names = F, sep = "\t",append=F)

  # subset the columns to only include the significant predictors (this removes the IDs)
  PRS_pheno_covs_cleaned = PRS_pheno_covs[,sig_predictors_cleaned]

  # if it contained the PRS, then we removed that too, as the PRS is passed in as a separate file if used at all
  if("SCORE1_SUM" %in% colnames(PRS_pheno_covs_cleaned)) {
    print("contains PRS, we remove it!")
    indexToDrop = which(colnames(PRS_pheno_covs_cleaned) == "SCORE1_SUM")
    PRS_pheno_covs_cleaned = PRS_pheno_covs_cleaned[,-indexToDrop ]
  }
  
  # split the covars into numeric and factor
  factorIndices = c()
  for(i in 1:length(colnames(PRS_pheno_covs_cleaned))) {
    
    if(grepl("fact", colnames(PRS_pheno_covs_cleaned)[i])) { # if its a factor, then it will have the "fact" in it
      factorIndices = c(factorIndices, i)
      print(paste0(colnames(PRS_pheno_covs_cleaned)[i], " is a factor"  ) )
    }
  }
  CovNumeric_toWriteOut = PRS_pheno_covs_cleaned[,-factorIndices]
  CovFactor_toWriteOut = PRS_pheno_covs_cleaned[,factorIndices]
  
  
  # write them out separately
  write.table(CovNumeric_toWriteOut,file=paste0(outL,"_","CovNumeric"), col.names = F, quote = F, row.names = F, sep = "\t",append=F)
  
  write.table(CovFactor_toWriteOut,file=paste0(outL,"_","CovFactor"), col.names = F, quote = F, row.names = F, sep = "\t",append=F)
  } # end of if we requested Covar writeout
}


# scenar = "PRS_Cov"
# scenar = "Cov"
# fit PRS + Cov
fitModel(train_data, test_data,"PRS_Cov")
dim(train_data)
scenar= "PRS_Cov"
print("now fitting covariates only model")
# refit model with just Cov by removing the PRS column
indexToDrop = which(colnames(train_data) == "SCORE1_SUM")
train_data = train_data[,-indexToDrop ]
test_data = test_data[,-indexToDrop ]

dim(train_data)
fitModel(train_data, test_data, "Cov", T)

#######################

# # refit model with only significant predictors
# train_data_reduced = train_data[, sig_predictors_cleaned]
# model2 = glm(train_pheno ~ .,family=binomial(link='logit'),data=train_data_reduced)
# test_data_reduced = test_data[, sig_predictors_cleaned]
# y_hat2 = predict(model2, newdata=test_data_reduced)
# cor(y_hat2, y_hat) # 0.974375, 0.9777882  (2nd val is when we fit lm instead of glm)
# cor(test_pheno, y_hat) #  0.3470513,  0.3525817
# cor(test_pheno, y_hat2) #  0.349472, 0.3502679
# # so they are very similar, and the reduced model appears to fit slightly better, but I think it is bad practice to refit model like that



###############################################
# train_fact_indices = c()
# colnam= colnames(train_data)
# for (i in 1:length(colnam)) {
#   if(grepl("fact", colnam[i])) { # if its a factor, then it will have the "fact" in it, as I've renamed all column names like that
#     #print( colnam[i])
#     train_fact_indices = c(train_fact_indices, i)
#   }
# }
# train_data_factors = train_data[,train_fact_indices]
# test_data_factors = test_data[,train_fact_indices]
# 
# indices_to_remove_train = sapply(train_data_factors, nlevels) <= 1 
# sum(indices_to_remove_train) #
# which(indices_to_remove_train == 1)
# # fact14_level
# # 11
# train_data_factors[,indices_to_remove_train]
# 
# table(train_data_factors$fact14_level)
# # 1
# # 31994
# table(test_data_factors$fact14_level)
# # 1     2
# # 10698     0
# table(PRS_pheno_covs$fact14_level)
# # 1     2
# # 53370     3
# # OK so I have introduced the redundant factor level, based on the valid/test split...
# # SOLUTION: I need to check it 3 ways, not just 2 ways...
# 
# 
# 
# indices_to_remove_test = sapply(test_data_factors, nlevels) <= 1 
# sum(indices_to_remove_test)
# which(indices_to_remove_test == 1)



# 
# 
# if( sum(indices_to_remove) > 0) {
#   colnamesToRemove = colnames(dataSplit)[indices_to_remove]
#   colindicesToRemove = match(colnamesToRemove, colnames(PRS_pheno_covs)) # as we subset the matrix, we need to map these back to the originals
#   PRS_pheno_covs = PRS_pheno_covs[,-colindicesToRemove]
# }



# CovFactors_used <- lapply(CovFactors_used, droplevels)
# CovFactors_used = as.data.frame(CovFactors_used)
# indices_to_remove = sapply(CovFactors_used, nlevels) <= 1    # remove factors with 1 level # https://stackoverflow.com/questions/17995195/remove-variables-with-factor-level-1
# 
# if( sum(indices_to_remove) > 0) {
#   colnamesToRemove = colnames(CovFactors_used)[indices_to_remove]
#   colindicesToRemove = match(colnamesToRemove, colnames(PRS_pheno_covs)) # as we subset the matrix, we need to map these back to the originals
#   PRS_pheno_covs = PRS_pheno_covs[,-colindicesToRemove]
# }
# 
# 
# # It is possible that we will create unused factor levels when splitting the data into training/test sets, in either of the training/test sets...
# # these levels need to be set to the mode in BOTH
# colnam= colnames(train_data)
# for (i in 1:length(colnam)) {
#   if(grepl("fact", colnam[i])) { # if its a factor, then it will have the "fact" in it, as I've renamed all column names like that
#     #print( colnam[i])
#     # tabulate levels
#     train_table = table(train_data[,colnam[i]])
#     test_table = table(test_data[,colnam[i]])
#     
#     # find any unused levels in either train/test
#     unused_train = as.numeric(which(train_table == 0))
#     unused_test = as.numeric(which(test_table == 0))
#     
#     # if there were any unused factors in the train, then we set those to be NA in the test
#     if ( length(unused_train) ) {
#       unusedFactors = unlist(dimnames(train_table))[unused_train]
#       medianvalue = names(sort(-table(test_data[,colnam[i]])))[1]#  median(test_data[,colnam[i]])
#       print(paste0("removing from Test data predictor ", colnam[i]," the unused train factors: ",unusedFactors, " and replacing it with mode: ", medianvalue ))
#       levels(test_data[,colnam[i]])[levels(test_data[,colnam[i]]) == unusedFactors] <- medianvalue # NA
#     }
#     
#     if ( length(unused_test) ) {
#       unusedFactors = unlist(dimnames(test_table))[unused_test]
#       medianvalue = names(sort(-table(train_data[,colnam[i]])))[1] #  median()
#       print(paste0("removing from Train data predictor ", colnam[i]," the unused test factors: ",unusedFactors, " and replacing it with mode: ", medianvalue ))
#       
#       levels(train_data[,colnam[i]])[levels(train_data[,colnam[i]]) == unusedFactors] <- medianvalue # NA
#     }
#   }
# }
# 



########################################
