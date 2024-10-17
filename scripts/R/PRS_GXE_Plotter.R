# creates final plots for the PRS GXE project

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")} 
 
baseLoc = args[1]
outputLoc= args[2]
plotName= args[3]

# Debug vars
# baseLoc =  "C:/0LocalHPC/PRS_GXE/240124/finalRes"
# outputLoc = "C:/0LocalHPC/PRS_GXE/240124/local/finalRes_res"
# plotName="PRS_GxE_non-linearity"

# baseLoc =  "C:/0LocalHPC/PRS_GXE/310124/finalRes"
# outputLoc = "C:/0LocalHPC/PRS_GXE/310124/local/finalRes_res"
# plotName="PRS_GxE_non-linearity"

################
# load the validation set performances too:
linear_valid= read.table(paste0(baseLoc,"_linear_valid"), header = T, sep="\t") 
nonlinear_valid= read.table(paste0(baseLoc,"_nonlinear_valid"), header = T, sep="\t") 


linear_valid_p= read.table(paste0(baseLoc,"_linear_valid"), header = T, sep="\t") 
nonlinear_valid_p= read.table(paste0(baseLoc,"_nonlinear_valid"), header = T, sep="\t") 

linear_valid_small= read.table(paste0(baseLoc,"_linear_valid_small"), header = T, sep="\t") 
nonlinear_valid_small= read.table(paste0(baseLoc,"_nonlinear_valid_small"), header = T, sep="\t") 

# get the indices for when the small is better for both the SNP.PRS and the SNP.PRS.Cov cases
nonlinear_smallBetter_SNPPRS = which (nonlinear_valid$SNP.PRS < nonlinear_valid_small$SNP.PRS) # 1  3  7 11 12 15 16 22 25
linear_smallBetter_SNPPRS = which (linear_valid$SNP.PRS < linear_valid_small$SNP.PRS)   
nonlinear_smallBetter_SNPPRSCOV = which (nonlinear_valid$SNP.PRS.Cov < nonlinear_valid_small$SNP.PRS.Cov) # 1  3  7 11 12 15 16 22 25
linear_smallBetter_SNPPRSCOV = which (linear_valid$SNP.PRS.Cov < linear_valid_small$SNP.PRS.Cov)   


# load in both the regular and the '_small' models, and keep the one that had better validation performance

# 1. load phenos and their respective p-values
linear= read.table(paste0(baseLoc,"_linear"), header = T, sep="\t") 
nonlinear= read.table(paste0(baseLoc,"_nonlinear"), header = T, sep="\t") 
linear_p= read.table(paste0(baseLoc,"_linear_p"), header = T, sep="\t") 
nonlinear_p= read.table(paste0(baseLoc,"_nonlinear_p"), header = T, sep="\t") 

linear_small= read.table(paste0(baseLoc,"_linear_small"), header = T, sep="\t") 
nonlinear_small= read.table(paste0(baseLoc,"_nonlinear_small"), header = T, sep="\t") 
linear_p_small= read.table(paste0(baseLoc,"_linear_p_small"), header = T, sep="\t") 
nonlinear_p_small= read.table(paste0(baseLoc,"_nonlinear_p_small"), header = T, sep="\t") 

# overwrite to use the better model's results
nonlinear$SNP.PRS[nonlinear_smallBetter_SNPPRS] = nonlinear_small$SNP.PRS[nonlinear_smallBetter_SNPPRS]
linear$SNP.PRS[linear_smallBetter_SNPPRS] = linear_small$SNP.PRS[linear_smallBetter_SNPPRS]
nonlinear_p$SNP.PRS[nonlinear_smallBetter_SNPPRS] = nonlinear_p_small$SNP.PRS[nonlinear_smallBetter_SNPPRS]
linear_p$SNP.PRS[linear_smallBetter_SNPPRS] = linear_p_small$SNP.PRS[linear_smallBetter_SNPPRS]

nonlinear$SNP.PRS.Cov[nonlinear_smallBetter_SNPPRSCOV] = nonlinear_small$SNP.PRS.Cov[nonlinear_smallBetter_SNPPRSCOV]
linear$SNP.PRS.Cov[linear_smallBetter_SNPPRSCOV] = linear_small$SNP.PRS.Cov[linear_smallBetter_SNPPRSCOV]
nonlinear_p$SNP.PRS.Cov[nonlinear_smallBetter_SNPPRSCOV] = nonlinear_p_small$SNP.PRS.Cov[nonlinear_smallBetter_SNPPRSCOV]
linear_p$SNP.PRS.Cov[linear_smallBetter_SNPPRSCOV] = linear_p_small$SNP.PRS.Cov[linear_smallBetter_SNPPRSCOV]

# Add an extra column to show if small model was used
nonlinear$smallModel_PRS.Cov = F
linear$smallModel_PRS.Cov = F
nonlinear$smallModel_PRS.Cov[nonlinear_smallBetter_SNPPRSCOV] =T
linear$smallModel_PRS.Cov[linear_smallBetter_SNPPRSCOV] = T

nonlinear$smallModel_SNP.PRS= F
linear$smallModel_SNP.PRS= F
nonlinear$smallModel_SNP.PRS[nonlinear_smallBetter_SNPPRS] = T
linear$smallModel_SNP.PRS[linear_smallBetter_SNPPRS] = T


# Standardise the above to be expressed as a %: 
nonlinear_perc = nonlinear[,c(2:ncol(nonlinear))]
linear_perc = linear[,c(2:ncol(linear))]

# 1. PRS results are expressed as units relative to the PRSIndi
nonlinear_perc$SNP.PRS = nonlinear$SNP.PRS / linear$PRSindi
linear_perc$SNP.PRS = linear$SNP.PRS / linear$PRSindi
nonlinear_perc$SNP = nonlinear$SNP / linear$PRSindi
linear_perc$SNP = linear$SNP / linear$PRSindi

# 2. PRS+Cov results are expressed as units relative to the PRSIndi+CovAdditive
#linear$PRSindi.CovAdd = abs(linear$PRSindi.CovAdd)
nonlinear_perc$SNP.PRS.Cov = nonlinear$SNP.PRS.Cov / linear$PRSindi.CovAdd
linear_perc$SNP.PRS.Cov = linear$SNP.PRS.Cov / linear$PRSindi.CovAdd


##################
# Find sensible subsets of the data, ie where the model has successfully converged and there was some signal
findNonOutliers_fixed = function(x, x_p, fixedThreshold = 5, alpha = 0.05) {
  nonOutliers = which(abs(x) < fixedThreshold & x_p <= alpha)
  return(nonOutliers)
}

findNonOutliers_baseline = function(x, x_p, alpha = 0.05) {
  meanx = mean(x) # 1.055619
  nonOutliers = which(x >  0 & x_p <= 0.05)
  x[nonOutliers]
  mean(x[nonOutliers])
  median(x[nonOutliers])
  return(nonOutliers)
}

linear_PRS = findNonOutliers_fixed(linear_perc$SNP.PRS, linear_p$SNP.PRS)
nonlinear_PRS = findNonOutliers_fixed(nonlinear_perc$SNP.PRS, nonlinear_p$SNP.PRS)
baseline_PRS = findNonOutliers_baseline(linear$PRSindi, linear_p$PRSindi)
PRS_allgood = intersect(intersect(linear_PRS,nonlinear_PRS),baseline_PRS)
#mean(linear_perc$SNP.PRS[PRS_allgood]) # 0.921359
#mean(nonlinear_perc$SNP.PRS[PRS_allgood]) # 0.9656543
#median(linear_perc$SNP.PRS[PRS_allgood]) #  0.8701804
#median(nonlinear_perc$SNP.PRS[PRS_allgood]) # 0.9319538

linear_SNP = findNonOutliers_fixed(linear_perc$SNP, linear_p$SNP)
nonlinear_SNP = findNonOutliers_fixed(nonlinear_perc$SNP, nonlinear_p$SNP)
baseline_SNP = findNonOutliers_baseline(linear$PRSindi, linear_p$PRSindi)
PRS_allgood_SNP = intersect(intersect(linear_SNP,nonlinear_SNP),baseline_SNP)
#mean(linear_perc$SNP[PRS_allgood_SNP]) # 0.7347101
#mean(nonlinear_perc$SNP[PRS_allgood_SNP]) #  0.631199
#median(linear_perc$SNP[PRS_allgood_SNP]) #  0.6155584
#median(nonlinear_perc$SNP[PRS_allgood_SNP]) # 0.5248266


linear_PRS.Cov = findNonOutliers_fixed(linear_perc$SNP.PRS.Cov, linear_p$SNP.PRS.Cov)
nonlinear_PRS.Cov = findNonOutliers_fixed(nonlinear_perc$SNP.PRS.Cov, nonlinear_p$SNP.PRS.Cov)
baseline_PRS.Cov = findNonOutliers_baseline(linear$PRSindi.Cov, linear_p$PRSindi.Cov)
PRS.Cov_allgood = intersect(intersect(linear_PRS.Cov,nonlinear_PRS.Cov),baseline_PRS.Cov)
#mean(linear_perc$SNP.PRS.Cov[PRS.Cov_allgood]) # 0.8563377
#mean(nonlinear_perc$SNP.PRS.Cov[PRS.Cov_allgood]) #  0.8943856
#median(linear_perc$SNP.PRS.Cov[PRS.Cov_allgood]) # 0.9137229
#median(nonlinear_perc$SNP.PRS.Cov[PRS.Cov_allgood]) # 0.9488377


# supp data 3 data:

# SNP.PRS
supp3_SNP.PRS = cbind.data.frame(linear$pheno, linear_perc$SNP.PRS , nonlinear_perc$SNP.PRS )
colnames(supp3_SNP.PRS) = c("pheno", "linear", "nonlinear")
supp3_SNP.PRS$linear = round(supp3_SNP.PRS$linear,3)
supp3_SNP.PRS$nonlinear = round(supp3_SNP.PRS$nonlinear,3)

# add stars where the small model was used
supp3_SNP.PRS$linear[linear_perc$smallModel_SNP.PRS] = paste0(supp3_SNP.PRS$linear[linear_perc$smallModel_SNP.PRS],"(s)")
supp3_SNP.PRS$nonlinear[nonlinear_perc$smallModel_SNP.PRS] = paste0(supp3_SNP.PRS$nonlinear[nonlinear_perc$smallModel_SNP.PRS],"(s)")




# SNP.Cov
supp3_SNP.Cov = cbind.data.frame(linear$pheno, linear_perc$SNP.PRS.Cov , nonlinear_perc$SNP.PRS.Cov )
colnames(supp3_SNP.Cov) = c("pheno", "linear", "nonlinear")
supp3_SNP.Cov$linear = round(supp3_SNP.Cov$linear,3)
supp3_SNP.Cov$nonlinear = round(supp3_SNP.Cov$nonlinear,3)

# add stars where the small model was used
supp3_SNP.Cov$linear[linear_perc$smallModel_PRS.Cov] = paste0(supp3_SNP.Cov$linear[linear_perc$smallModel_PRS.Cov],"(s)")
supp3_SNP.Cov$nonlinear[nonlinear_perc$smallModel_PRS.Cov] = paste0(supp3_SNP.Cov$nonlinear[nonlinear_perc$smallModel_PRS.Cov],"(s)")

# subset to actually used data:
supp3_SNP.PRS = supp3_SNP.PRS[PRS_allgood,]
supp3_SNP.Cov = supp3_SNP.Cov[PRS.Cov_allgood,]

write.table(supp3_SNP.PRS,paste0(outputLoc,"_supp3_SNPPRS.txt" ), quote = F, row.names = F, sep = "\t" )
write.table(supp3_SNP.Cov,paste0(outputLoc,"_supp3_SNPCovS.txt" ), quote = F, row.names = F, sep = "\t" )
nrow(supp3_SNP.PRS)
nrow(supp3_SNP.Cov)
#median(supp3_SNP.PRS$linear) # 0.87
#median(supp3_SNP.PRS$nonlinear) # 0.932
#median(supp3_SNP.Cov$linear) # 0.914
#median(supp3_SNP.Cov$nonlinear) #  0.949



# Final data:
NN.SNP_linear = linear_perc$SNP[PRS_allgood_SNP]
NN.SNP_nonlinear = nonlinear_perc$SNP[PRS_allgood_SNP]
NN.SNP_linear_raw = linear$SNP[PRS_allgood_SNP]
NN.SNP_nonlinear_raw = nonlinear$SNP[PRS_allgood_SNP]
#t.test(NN.SNP_linear_raw,NN.SNP_nonlinear_raw, paired =T)$p.value # 0.5293173

NN.SNP.PRS_linear = linear_perc$SNP.PRS[PRS_allgood]
NN.SNP.PRS_nonlinear = nonlinear_perc$SNP.PRS[PRS_allgood]
NN.SNP.PRS_linear_raw = linear$SNP.PRS[PRS_allgood]
NN.SNP.PRS_nonlinear_raw = nonlinear$SNP.PRS[PRS_allgood]
#t.test(NN.SNP.PRS_linear_raw,NN.SNP.PRS_nonlinear_raw, paired =T)$p.value # 0.07809359


NN.SNP.PRS.Cov_linear = linear_perc$SNP.PRS.Cov[PRS.Cov_allgood]
NN.SNP.PRS.Cov_nonlinear = nonlinear_perc$SNP.PRS.Cov[PRS.Cov_allgood]
#length(NN.SNP.PRS.Cov_linear)
NN.SNP.PRS.Cov_linear_raw = linear$SNP.PRS.Cov[PRS.Cov_allgood]
NN.SNP.PRS.Cov_nonlinear_raw = nonlinear$SNP.PRS.Cov[PRS.Cov_allgood]
#t.test(NN.SNP.PRS.Cov_linear_raw,NN.SNP.PRS.Cov_nonlinear_raw, paired =T)$p.value # 0.8762695

#t.test( c(NN.SNP_linear_raw, NN.SNP.PRS_linear_raw, NN.SNP.PRS.Cov_linear_raw),c(NN.SNP_nonlinear_raw, NN.SNP.PRS_nonlinear_raw, NN.SNP.PRS.Cov_nonlinear_raw), paired =T)$p.value # 0.5827469



# final table of medians
finalTable = cbind(median (NN.SNP_linear), median(NN.SNP.PRS_linear), median(NN.SNP.PRS.Cov_linear) )
finalTable= rbind(finalTable, cbind(median(NN.SNP_nonlinear), median(NN.SNP.PRS_nonlinear), median(NN.SNP.PRS.Cov_nonlinear)))
finalTable= rbind(finalTable, cbind(t.test(NN.SNP_linear,NN.SNP_nonlinear, paired = T)$p.value, t.test(NN.SNP.PRS_linear,NN.SNP.PRS_nonlinear, paired = T)$p.value, t.test(NN.SNP.PRS.Cov_linear,NN.SNP.PRS.Cov_nonlinear, paired = T)$p.value))

colnames(finalTable) = c("SNP", "SNP+PRS", "SNP+PRS+Cov")
rownames(finalTable) = c("linear", "non-linear", "t-test p")
finalTable = finalTable[1:2,2:3]
finalTable = as.data.frame(finalTable)

boxLabels = c(finalTable$`SNP+PRS`[1], finalTable$`SNP+PRS`[2], finalTable$`SNP+PRS+Cov`[1], finalTable$`SNP+PRS+Cov`[2])
boxLabels = round(boxLabels,3)
write.table(round(finalTable,3),paste0(outputLoc,".txt" ), quote = F )

# final table of means
finalTable = cbind(mean (NN.SNP_linear), mean(NN.SNP.PRS_linear), mean(NN.SNP.PRS.Cov_linear) )
finalTable= rbind(finalTable, cbind(mean(NN.SNP_nonlinear), mean(NN.SNP.PRS_nonlinear), mean(NN.SNP.PRS.Cov_nonlinear)))
finalTable= rbind(finalTable, cbind(t.test(NN.SNP_linear,NN.SNP_nonlinear, paired = T)$p.value, t.test(NN.SNP.PRS_linear,NN.SNP.PRS_nonlinear, paired = T)$p.value, t.test(NN.SNP.PRS.Cov_linear,NN.SNP.PRS.Cov_nonlinear, paired = T)$p.value))

colnames(finalTable) = c("SNP", "SNP+PRS", "SNP+PRS+Cov")
rownames(finalTable) = c("linear", "non-linear", "t-test p")
finalTable = finalTable[1:2,2:3]
finalTable = as.data.frame(finalTable)


write.table(round(finalTable,3),paste0(outputLoc,"_mean.txt" ), quote = F )


####################

# plot main results
colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
allColours=colorBlind7[1:2] # pick 4 nice colors for the 4 estimates (Vg, Ve, VGxE, Fenv)

plotName= "" # plotname disabled for publication

# find the nearest 0.5 of all of the data
alldata = c(NN.SNP.PRS_linear, NN.SNP.PRS_nonlinear, NN.SNP.PRS.Cov_linear, NN.SNP.PRS.Cov_nonlinear)
minVal = round(min(alldata)*2) / 2 #  = -0.5
maxVal = round(max(alldata)*2) / 2 # + 0.25 # = 4  + 0.2 for the text labels

boxNames = c("SNP", "SNP+Cov") # each box will have the name of the model (ie what we are testing, NOT the truth)
pixelScale = 100
yaxislab = "fraction of baseline"


allMeans = list(median(NN.SNP.PRS_linear), median(NN.SNP.PRS_nonlinear), median(NN.SNP.PRS.Cov_linear), median(NN.SNP.PRS.Cov_nonlinear))
inputDatas_mat = list(NN.SNP.PRS_linear, NN.SNP.PRS_nonlinear, NN.SNP.PRS.Cov_linear, NN.SNP.PRS.Cov_nonlinear)

# insert a blank column in the middle to separate out the 2 scenarios: https://stackoverflow.com/questions/58708529/how-to-create-spaces-between-groups-and-control-size-of-axis-labels-in-boxplot
VEC = 1:4
VEC[3:length(VEC)] = VEC[3:length(VEC)] +1

filen = paste(outputLoc,".png", sep="" )

png(filen, width=3.5* length(boxNames) * pixelScale , height=4.4 * pixelScale);
xvals = c(1,2,4,5)
xvals = xvals + 0.3 # so that it is slightly offset to the right so it is visible
labelCols  = c(allColours,allColours)
i=1

boxplot(allMeans , outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0,xaxt="n", border =allColours, xlab="", ylab=yaxislab, main=plotName, ylim=c(minVal, maxVal), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8, at=VEC)
stripchart( inputDatas_mat, vertical = TRUE, pch = 16, col = allColours, add = TRUE, method = "jitter", cex = 2,at=VEC) #

# add the median values at the top of each stripchart
for(i in 1:length(boxLabels)) {
  yval = allMeans[[i]] + 0.25# max(inputDatas_mat[[i]]) + 0.25
  text(xvals[i], yval, boxLabels[i], cex = 0.88, col = labelCols[i])
 # print(paste0("adding label ", boxLabels[i], " to x: ", xvals[i], " y: ", yval ))
}


axis(1,at=c(1.5, 4.5),labels=boxNames, cex.axis = 2,mgp=c(3,1.5,0)) 
legend("top", legend=c("linear NN","nonlinear NN"),col=allColours, text.col= allColours, lty=1, lwd=3, cex = 1.5)
title(xlab = "scenario", line = 3.7, cex.lab = 1.5)

dev.off()


# calculate % differences too
percDiff = function (a, b, dec = 4) {
  return( round(  abs(a - b) / ((a + b) / 2) , dec) )
}


# median
print (paste0 ("SNP linear to non-linear % difference: ", percDiff(median(NN.SNP.PRS_linear), median(NN.SNP.PRS_nonlinear)),  " of ", length(NN.SNP.PRS_linear) , " traits"   ) )
# "SNP linear to non-linear % difference: 0.0686 of 27 traits"



print (paste0 ("SNP+Cov linear to non-linear % difference: ", percDiff(median(NN.SNP.PRS.Cov_linear), median(NN.SNP.PRS.Cov_nonlinear)),  " of ", length(NN.SNP.PRS.Cov_nonlinear) , " traits"   )  ) 
# "SNP+Cov linear to non-linear % difference: 0.0377 of 21 traits"


print (paste0 ("SNP non-linear to baseline % difference: ", percDiff(1, median(NN.SNP.PRS_nonlinear))  ) )
# "non-linear to baseline % difference: 0.0704"

print (paste0 ("SNP+Cov non-linear to baseline % difference: ", percDiff(1, median(NN.SNP.PRS.Cov_nonlinear))  ) )
# "SNP+Cov non-linear to baselin e% difference: 0.0525"



# So it seems that there is again, some, minimal amount of non-linearity, but it is NOT better than the basic, additive regression model

###############################################
# Print out the median fractions too
# (these will be a few decimals different from the 1-above. This is because above we use the percDiff, whereas now we just use 1-fraction)

print (paste0 ("SNP non-linear fraction: ", round(median(NN.SNP.PRS_nonlinear), 4),  " of ", length(NN.SNP.PRS_linear) , " traits"   ) )
# "SNP non-linear fraction: 0.932 of 27 traits"
print (paste0 ("SNP linear fraction: ",round(  median(NN.SNP.PRS_linear), 4),  " of ", length(NN.SNP.PRS_linear) , " traits"   ) )
# "SNP linear fraction: 0.8702 of 27 traits"


print (paste0 ("SNP+Cov non-linear fraction: ", round( median(NN.SNP.PRS.Cov_nonlinear), 4),  " of ", length(NN.SNP.PRS.Cov_nonlinear) , " traits"   )  ) 
# "SNP+Cov non-linear fraction: 0.9488 of 21 traits"

print (paste0 ("SNP+Cov linear fraction: ", round(  median(NN.SNP.PRS.Cov_linear), 4),  " of ", length(NN.SNP.PRS.Cov_nonlinear) , " traits"   )  ) 
# "SNP+Cov linear fraction: 0.9137 of 21 traits"

#round(0.931953831026769, 4)






