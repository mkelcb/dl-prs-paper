# correlates the 2nd and 3rd cols of a csv and append writes the r^2 to file
# also optionally calculates AUC

library(pROC)
library(r2redux)
library(largesamplehl)

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 
 
gwasLoc = args[1]
outputLoc= args[2]
calcAUC = args[3] == "1"
hasHeader = is.na(args[4]) == F && args[4] == "1"
doAppend = is.na(args[5]) == F && args[5] == "1"
separator= ","
if ( is.na(args[6]) == F && args[6] == "1") { separator = "\t"}

print(paste0("hasHeader: ", hasHeader, " / doAppend: " , doAppend , " / calcAUC: ", calcAUC, " / separator: ", separator))

# Debug vars
# gwasLoc =  "C:/0Datasets/shaPRS/crossAncestry/RapidoPGS_EUR_height_raw_noAmbiguousAlleles.csv" # "C:/0Datasets/shaPRS/crossAncestry/RapidoPGS_JP_height_raw_noAmbiguousAlleles.csv"  "C:/0Datasets/shaPRS/crossAncestry/RapidoPGS_EUR_JAP_height_max_combined.csv"  "C:/0Datasets/shaPRS/crossAncestry/RapidoPGS_EUR_JAP_height_max.csv"  "C:/0LocalHPC/results/AAD_any_all.sscore"
# outputLoc =   "C:/0Datasets/shaPRS/crossAncestry/RapidoPGS_EUR_JAP_height_max"
#hasHeader= FALSE
#doAppend= FALSE 
#calcAUC=FALSE 
#separator= ","


# 1. load phenos
gwasRes= read.table(gwasLoc, header = hasHeader, sep=separator) 
gwasRes =na.omit(gwasRes)

# also remove -9s as those indicate missing pheno for a .fam file
gwasRes = gwasRes[gwasRes[,2] != -9, ]


correlation = cor(gwasRes[,2], gwasRes[,3])
correlation_sq = correlation^2
test = cor.test(gwasRes[,2], gwasRes[,3])

# get confidence intervals via r2redux
dat1 =  cbind(gwasRes[,2], gwasRes[,3])
nv= nrow(dat1)
output1=r2_var(dat1,1,nv)
#?r2_var
#output1$r2_based_p # why this is DIFFERENT? ??
#test$p.value

testResults = cbind.data.frame(correlation,correlation_sq,test$p.value, output1$r2_based_p, output1$lower_r2, output1$upper_r2	)
colnames(testResults) = c("correlation","correlation^2","corr_p","r2redux_p", "CI_LOW", "CI_HIGH")
if(doAppend == F) {
  write.table(testResults,file=outputLoc, quote = F, row.names = F, sep = "\t",append=doAppend)
} else { # if we are appending only write the squared correlations, as we assume that will be a column of data
  write.table(correlation_sq,file=outputLoc, quote = F,col.names = F, row.names = F, sep = "\t",append=doAppend)
  write.table(testResults,file=paste0(outputLoc,"_redux"), quote = F,col.names = F, row.names = F, sep = "\t",append=doAppend)
  
}


print(paste("written: correlation_sq",signif(correlation_sq,3), "(CI_LOW:",signif(output1$lower_r,3), "CI_HIGH:",signif(output1$upper_r2,3),") r2_based_p: ", signif(output1$r2_based_p,3)   , "to", outputLoc ))


if(calcAUC) {
  
  ###########
  # Perform modified HL test:  https://github.com/mbounthavong/R-tutorials/blob/main/Logistic%20Regression%20in%20R-Part2
  binaryOutcome= gwasRes[,2] -1 # outcome should be between 0 and 1
  PRS= gwasRes[,3]
  
  #Fit a logistic model, regressing the true pheno on the PRS, https://github.com/nikola-sur/goodGLM
  model <- glm(binaryOutcome ~ PRS, family = binomial(link="logit"))
  hltest_rest = hltest(model)
  
  
  # same as above
  # phat <- predict(model, type = "response")
  # hltest(y=binaryOutcome, prob=phat) 
  
  ##############
  
  

category= gwasRes[,2]
prediction= gwasRes[,3]
# category <- c(1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
# prediction <- rev(seq_along(category))
#prediction[9:10] <- mean(prediction[9:10])
roc_obj <- roc(category, prediction)
myAUC= auc(roc_obj)[1]

ci_res <- ci(roc_obj)

#data(aSAH)
#roc_obj <- roc(aSAH$outcome, aSAH$s100b)


# AUC
print("cases / control table")
as.data.frame(table(gwasRes[,2]))
filen= paste0(outputLoc, "_AUC")
AUCTable=cbind(myAUC[1],ci_res[1] ,ci_res[3], hltest_rest$p.value)
colnames(AUCTable) = c("AUC", "CI_LOW", "CI_HIGH", "HL_p")



if(doAppend == F) {
  write.table(AUCTable,file=filen, quote = F, row.names = F, sep = "\t")
} else { # if we are appending only write the squared correlations, as we assume that will be a column of data
  write.table(myAUC[1],file=filen, quote = F,col.names = F, row.names = F, sep = "\t",append=doAppend)
  
}




print(paste("written: AUC",signif(myAUC,3), "CI low:",signif(ci_res[1],4), " / high: ",signif(ci_res[3],4), "/ modified HL test p:",  signif(hltest_rest$p.value,4)  , "to", filen ))

filename=paste0(outputLoc, "_AUC.png")
png(filename, width=640 , height=640, res =128);

#plot.roc(category,prediction, ci=TRUE, of="thresholds")
plot.roc(category,prediction)
dev.off()



}


####################
