# creates 'endpoints' directory structure from the PGSlist_50K_subset_annotated.csv


args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")} 


dataLoc= args[1]
outputLoc= args[2]

#DEBUG VARS
#dataLoc= "C:/Users/M Kel/GoogleDrive_Cam/0Publications/PRS_GXE/data/PGSlist_50K_subset_annotated.csv"
#outputLoc= "C:/Users/M Kel/GoogleDrive_Cam/0Publications/PRS_GXE/data/endpoints/"

phenoTextList = "phearray=( "
phenoTextList_all = "phearray_withCEUdef=( "
binaryTextList = "phearray_binary=( "
sex_subsetList = "sex_subsetList=( "
covariateExclusionList = "covariateExclusionList=( "
data  = read.table(dataLoc, header = T, sep="\t", quote = "\"")
for (i in 1:nrow(data)) {
  print(paste0("processing pheno: ", data$trait[i]))
  phenoTextList_all = paste0(phenoTextList_all,"'" ,data$EFO[i],"'" ," ")
  sex_subsetList = paste0(sex_subsetList,"'" ,data$Sex_subset[i],"'" ," ")
  covariateExclusionList = paste0(covariateExclusionList,"'" ,data$covariate_exclusion[i],"'" ," ")
  binaryTextList = paste0(binaryTextList,"'" ,data$isBinary[i],"'" ," ")
  
  
  if(data$isEndpoint[i] == T) {
  phenoTextList = paste0(phenoTextList,"'" ,data$EFO[i],"'" ," ")
  txt = paste0("# ",data$trait[i], " definition:\n")
  if(data$ICD10[i] != "") {txt = paste0(txt, "ICD-10: ",data$ICD10[i], "\n")}
  if(data$ICD9[i] != "") {txt = paste0(txt, "prevalent icd-9: ",data$ICD9[i], "\n")}
  if(data$OPCS4[i] != "") {txt = paste0(txt, "OPCS-4: ",data$OPCS4[i], "\n")}
  if(data$OPCS3[i] != "") {txt = paste0(txt, "OPCS-3: ",data$OPCS3[i], "\n")}
  if(data$field.id[i] != "") {txt = paste0(txt, "Self-report field ",data$field.id[i],": ",data$code[i], "\n")}
  dir.create(file.path(paste0(outputLoc, data$EFO[i])) )
  #fileConn<-file(paste0(outputLoc, data$EFO[i],"/endpoint_definition.txt"))
  #writeLines(txt, fileConn)
  #close(fileConn)
  }
  
}

# print out the full list of phenos and other files, which are always aligned to the data
# these are to be pasted into the PRS_GX_functions_v2.sh and loaded as a session variable
phenoTextList = paste0(phenoTextList, ")")
print(phenoTextList)


phenoTextList_all = paste0(phenoTextList_all, ")")
print(phenoTextList_all)

binaryTextList = paste0(binaryTextList, ")")
print(binaryTextList)

sex_subsetList = paste0(sex_subsetList, ")")
print(sex_subsetList)


covariateExclusionList = paste0(covariateExclusionList, ")")
print(covariateExclusionList)




