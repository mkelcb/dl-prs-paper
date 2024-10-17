# creates final plots for the Simulations of the PRS GXE project

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")} 
 
baseLoc = args[1]
outputLoc= args[2]
plotName= args[3]

# Debug vars
# baseLoc =  "C:/0LocalHPC/PRS_GXE/simRes"
# outputLoc = "C:/0LocalHPC/PRS_GXE/simRes_res"
# plotName="PRS_GxE_simulations"

linear= read.table(paste0(baseLoc,"_linear"), header = T, sep="\t", fill = TRUE ) 



# 1. load phenos
linear= read.table(paste0(baseLoc,"_linear"), header = T, sep="\t") 
nonlinear= read.table(paste0(baseLoc,"_nonlinear"), header = T, sep="\t") 
noAct= read.table(paste0(baseLoc,"_noAct"), header = T, sep="\t") 

# remove 1st col, as that is the phenotype/replica counter
linear_perc = linear[,c(2:ncol(linear))]
nonlinear_perc = nonlinear[,c(2:ncol(nonlinear))] 
noAct_perc = noAct[,c(2:ncol(noAct))] 

# also get difference BEFORE standardising, so that we could use this to get the p.vals on the original scale
nonlinear_linear_diff_orig = nonlinear_perc - linear_perc

# standardise it to be a % of the additive PRS
#linear_perc = linear_perc / linear_perc$PRSindi
#nonlinear_perc = nonlinear_perc / nonlinear_perc$PRSindi
#noAct_perc = noAct_perc / noAct_perc$PRSindi
# dont do the above, we want to standardise each PRS with the correct baseline

# standardise against the additive phenotype's regular PRS
linear_perc$res = linear_perc$res / linear_perc$PRSIndi_1
nonlinear_perc$res = nonlinear_perc$res / nonlinear_perc$PRSIndi_1
noAct_perc$res = noAct_perc$res / noAct_perc$PRSIndi_1
linear_perc$res_PRS = linear_perc$res_PRS / linear_perc$PRSIndi_1
nonlinear_perc$res_PRS = nonlinear_perc$res_PRS / nonlinear_perc$PRSIndi_1
noAct_perc$res_PRS = noAct_perc$res_PRS / noAct_perc$PRSIndi_1

# standardise against the additive phenotype's haplotype effect PRS
linear_perc$res_haplo = linear_perc$res_haplo / linear_perc$PRSIndi_1_haplo
nonlinear_perc$res_haplo = nonlinear_perc$res_haplo / nonlinear_perc$PRSIndi_1_haplo
noAct_perc$res_haplo = noAct_perc$res_haplo / noAct_perc$PRSIndi_1_haplo
linear_perc$res_haplo_PRS = linear_perc$res_haplo_PRS / linear_perc$PRSIndi_1_haplo
nonlinear_perc$res_haplo_PRS = nonlinear_perc$res_haplo_PRS / nonlinear_perc$PRSIndi_1_haplo
noAct_perc$res_haplo_PRS = noAct_perc$res_haplo_PRS / noAct_perc$PRSIndi_1_haplo
linear_perc$res_haplo_filtered = linear_perc$res_haplo_filtered / linear_perc$PRSIndi_1_haplo
nonlinear_perc$res_haplo_filtered = nonlinear_perc$res_haplo_filtered / nonlinear_perc$PRSIndi_1_haplo
noAct_perc$res_haplo_filtered = noAct_perc$res_haplo_filtered / noAct_perc$PRSIndi_1_haplo

# standardise against the 4-way epistasis' phenotype's regular PRS
linear_perc$res_epi = linear_perc$res_epi / linear_perc$PRSIndi_4
nonlinear_perc$res_epi = nonlinear_perc$res_epi / nonlinear_perc$PRSIndi_4
noAct_perc$res_epi = noAct_perc$res_epi / noAct_perc$PRSIndi_4
linear_perc$res_epi_PRS = linear_perc$res_epi_PRS / linear_perc$PRSIndi_4
nonlinear_perc$res_epi_PRS = nonlinear_perc$res_epi_PRS / nonlinear_perc$PRSIndi_4
noAct_perc$res_epi_PRS = noAct_perc$res_epi_PRS / noAct_perc$PRSIndi_4




nonlinear_linear_diff = nonlinear_perc - linear_perc

#colnames(nonlinear_linear_diff) = c("res_epi_nonlinear", "res_nonlinear", "res_haplo_nonlinear", "res_haplo_PRS_nonlinear", "res_haplo_filtered_nonlinear", "a", "b")




# Plot the following 5 columns
# nonlinear_linear_diff$res_epi
# nonlinear_linear_diff$res
# nonlinear_linear_diff$res_haplo
# nonlinear_linear_diff$res_haplo_PRS
# nonlinear_linear_diff$res_haplo_filtered
plotData = nonlinear_linear_diff[,c("res_epi", "res", "res_haplo", "res_haplo_PRS", "res_haplo_filtered", "res_epi_PRS")]
plotData_orig = nonlinear_linear_diff_orig[,c("res_epi", "res", "res_haplo", "res_haplo_PRS", "res_haplo_filtered", "res_epi_PRS")]

# get a paired t-test result for each  ( a paired t-test is equivalent to the difference between the 2 vectors, so we are already good)
tableData = NULL
for (i in 1:ncol(plotData) ) {
  currData = plotData[,i]
  currData_orig = plotData_orig[,i]
  res2 <- t.test(currData) # we should use the ORIGINAL table for the t.test
  res2 <- t.test(currData_orig) # we should use the ORIGINAL table for the t.test
  
  tableData = rbind(tableData, cbind( colnames(plotData)[i] , mean(currData), res2$p.value ) )
  
}

cols= c("scenario", "mean_nonlinear", "t-test_p")
colnames(tableData) = cols
write.table(tableData,file=paste0(outputLoc,"_tab"), quote = F, row.names = F, sep = "\t",append=F)


res_PRS_adv = nonlinear_perc$res_PRS - nonlinear_perc$res	
res_PRS_adv_p <- t.test(res_PRS_adv)
res_epi_PRS_adv = nonlinear_perc$res_epi_PRS - nonlinear_perc$res_epi	
res_epi_PRS_adv_p <- t.test(res_epi_PRS_adv)
myDat = cbind("res_PRS_adv", mean(res_PRS_adv), res_PRS_adv_p$p.value)
myDat = rbind(myDat, cbind("res_epi_PRS_adv",mean(res_epi_PRS_adv), res_epi_PRS_adv_p$p.value) )
colnames(myDat) = cols
write.table(myDat,file=paste0(outputLoc,"_adv_tab"), quote = F, row.names = F, sep = "\t",append=F)

#######################################

dotPlotDataFrame= function (nonlinear_linear_diff, outputLoc) {
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77","#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888") # https://emitanaka.org/blog/2022-02-20-color-considerations/color-considerations.html#fig:safe
  all_box_names = c("a", "b", "c", "d", "e","f","g","h", "i", "j", "","k")
  
  allColours = safe_colorblind_palette[1:ncol(nonlinear_linear_diff)] # get some nice distinct colours based on the number of elements we want to plot
  allColours[allColours== "#FFFF00"] <- "black" # yellow is poor for reading with letters, replace it

  i=1
  # load ^ process data
  allData = NULL
  allMeans= list() # = vector(length = length(inputFiles) )
  boxNames = vector(length = ncol(nonlinear_linear_diff) )
  for (i in 1:ncol(nonlinear_linear_diff)) {
    data = nonlinear_linear_diff[,i]
    numNAs = length( which(is.na(data)) )
    if(numNAs > 0 ) { print("Input has NAs which were replaced by column mean")}
    data[is.na(data)] <- mean(data, na.rm = TRUE) # replace NA's by data mean
    
    allData = cbind(allData, data )
    allMeans[[i]] = mean( data )
    boxNames[i] = colnames(nonlinear_linear_diff)[i]
  }
  # actually just use alphabet for boxnames..
  boxNames = all_box_names[1:ncol(nonlinear_linear_diff)]
  
  allData[is.na(allData)] <- 0

  overallMean = mean( unlist(allMeans) )
  inputDatas_mat = list()
  for (i in 1:ncol(nonlinear_linear_diff)) {
    inputDatas_mat[[i]] = unlist(nonlinear_linear_diff[,i]) 
  }
  
  
  
  plotName ="" #  gsub("_", " ", plotName) #remove underscores
  plotNames= colnames(nonlinear_linear_diff)
  
  k=2
  u=1
  for (u in 1: ncol(nonlinear_linear_diff)) {
    minVal = min(inputDatas_mat[[u]]) * 0.9
    maxVal = max(inputDatas_mat[[u]]) * 1.1
    
    for (k in 1:2) {
      #plotName="" # plotname disabled for publication
      if( k == 1) {
      filen =paste0(outputLoc,plotNames[u],".pdf" )
      pdf(filen, width=3.8 , height=6.4);
      } else {
        filen =paste0(outputLoc,plotNames[u],".png" )
        png(filen, width=258.75  , height=200*3);
      }
      
      
      #par( mfrow = c( 1, 1 ) ,mgp = c(1, 2.5, 0))
      borderEnabled=allColours[u]
      op <- par(mar = c(5,5,4,2) + 0.1)
      boxplot(allMeans[u] , outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0,xaxt="n", border =borderEnabled, names=boxNames[u], ylab=  expression(paste("r"^"2")) , main=plotName, cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8, ylim=c(minVal, maxVal)) # 
      stripchart( inputDatas_mat[[u]], vertical = TRUE, pch = 16, col = allColours[u], add = TRUE, method = "jitter", cex = 2) #
      axis(1,at=c(1:1),labels=boxNames[u], cex.axis = 2,mgp=c(3,3.5,0))
      
      
      dev.off()
    }
  
  }
}
################


dotPlotDataFrame(nonlinear_linear_diff,outputLoc)


###############

# before <- linear$SNP.PRS.Cov
# after <- nonlinear$SNP.PRS.Cov
# res <- t.test(before, after, paired = TRUE)
# res$p.value
# 
# before_after = linear$SNP.PRS.Cov - nonlinear$SNP.PRS.Cov
# res2 <- t.test(before_after)
# res2$p.value
