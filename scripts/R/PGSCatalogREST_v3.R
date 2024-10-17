# https://www.geeksforgeeks.org/accessing-rest-api-using-r-programming/

# Installing the packages
#install.packages("httr")
#install.packages("jsonlite")

# Loading packages
library(httr)
library(jsonlite)

# Initializing API Call
numPGS= 3688 # but there are some redundant ones:

# https://www.pgscatalog.org/rest/score/PGS003764
# this is the maximum, by chcking it on the website
# also need to filter out ones that do NOT have the cohorts specified at all, as those could be anything
outLoc="C:/Users/M Kel/GoogleDrive_Cam/0Publications/PRS_GXE/data/PGSlist_50K"
currentNum = 1
maxPGSnum = 9999
minSampleSize = 50000  # 100000
allResults = NULL


#currentNum = PGS002978

#for (currentNum in 1:maxPGSnum) {
currentNum=1
while (currentNum < maxPGSnum) {
print(paste0("processing PGS num: ", currentNum))
PGS00 = paste0("PGS00", sprintf("%04d", currentNum))
call <- paste0("https://www.pgscatalog.org/rest/score/",  PGS00 )  #"https://www.pgscatalog.org/rest/score/PGS017889" # "https://www.pgscatalog.org/rest/score/PGS000013"
currentNum = currentNum +1


# Getting details in API
get_pgs_details <- GET(url = call)

# Getting status of HTTP Call
status_code(get_pgs_details) # this is 200, even when there is no PGS

# Content in the API
content = content(get_pgs_details)

# scraping a web service is stochastic, we may not get a result every time, even if there is data, so we try again 50 times
if(is.null(content$message) == F && content$message == "request limit exceeded") {
  pause = as.numeric( unlist( strsplit(content$availableIn," ") )[1] )
  
  print(paste0("request limit exceeded, retry in: ", pause))
  Sys.sleep(pause + 5) # add 5 secs to be sure
  
  
  #if( tryAgainCount <tryAgainMax) {
  #  print(paste0("No content for PGS num: ", currentNum, " trying again: ",tryAgainCount ))
  #  tryAgainCount = tryAgainCount+1
    currentNum = currentNum-1 # we go back one, to try the same again
  #} else {tryAgainCount = 0 } # edge case if we run out of 'retries' we must still reset the tryagain count
  next
}
tryAgainCount = 0 # reset this


if ( length(content) == 0 ) {
  

  if(currentNum > numPGS) {break} else {next}   
  
  } # if there is no content, that means we run out of PGS, we break, or just continue, if we still havent found most (as there are some 'gaps')


# Converting content to text
get_pgs_text <- content(get_pgs_details,"text", encoding = "UTF-8")

# Parsing data in JSON
get_pgs_json <- fromJSON(get_pgs_text,flatten = TRUE)
EFO = get_pgs_json$trait_efo$id




sampletypes = c('samples_variants', 'samples_training')
qcohort = 'UKB'
PGSNotEligible = F
numMissingCohortInfos = 0
for (sampletype in sampletypes) {
  for (sample in get_pgs_json[sampletype]){
    # Check cohorts
    if (length(sample) > 0) {
      for (cohort in sample['cohorts']){
        cohort = unlist(cohort)
        if ( length(cohort) > 0 ) {
        #print(paste0("cohort[1] is: ", cohort[1]))
        #if (cohort['name_short'] == qcohort) { PGSNotEligible = T   } # this doesnt work for some reason
        if (cohort[1] == qcohort) { PGSNotEligible = T   }
        } else {numMissingCohortInfos = numMissingCohortInfos +1} # if there was no data then, we take note
      }
    } else {numMissingCohortInfos = numMissingCohortInfos +1} # if there was no data then, we take note

  }
}

# actually, also exclude PGS that have interaction terms
if( is.null( get_pgs_json$variants_interactions) || get_pgs_json$variants_interactions > 0 )  {PGSNotEligible = T}

ancestryIsEuropean = F
if (length(get_pgs_json$samples_variants) > 0 && get_pgs_json$samples_variants$ancestry_broad == "European") {ancestryIsEuropean = T}
if (length(get_pgs_json$samples_training) > 0 && get_pgs_json$samples_training$ancestry_broad == "European") {ancestryIsEuropean = T}

sampleSize = 0
if (length(get_pgs_json$samples_variants) > 0 ) {sampleSize =                 sum(get_pgs_json$samples_variants$sample_number) }
if (length(get_pgs_json$samples_training) > 0 ) {sampleSize = max(sampleSize, sum(get_pgs_json$samples_training$sample_number)) }

PMID = get_pgs_json$publication$PMID
link = get_pgs_json$ftp_harmonized_scoring_files$GRCh37$positions
# if it included the UKBB, OR it had no information on both slots, then we don't trust it, it also must be European, and not yet have it yet
if ( numMissingCohortInfos < 2 && PGSNotEligible == F && ancestryIsEuropean && sampleSize >= minSampleSize && link %in% allResults[,4] == F) {
  
  
  trait = get_pgs_json$trait_reported
  
  
  
  method_name = get_pgs_json$method_name
  variants_number = get_pgs_json$variants_number
  
  newEntry = cbind(trait, sampleSize, PMID,link, method_name, variants_number, PGS00, EFO)
  allResults = rbind(allResults, newEntry)
  print(paste0("PGS", currentNum," (",trait,") met criteria"))

}

}

#allResults2 = read.table(paste0(outLoc, "_ALL.txt"), header = T, sep = "\t")
write.table(allResults, paste0(outLoc, "_ALL.txt"), sep = "\t", row.names = F, col.names = T, quote = T) # must use quote, as various characters break the columns


# manually unify spelling variations
# allResults[which(allResults[,1] == "Type II diabetes"),1] = "Type 2 diabetes"
# allResults[which(allResults[,1] == "Type 2 Diabetes"),1] = "Type 2 diabetes"
# allResults[which(allResults[,1] == "type 2 diabetes"),1] = "Type 2 diabetes"
# allResults[which(allResults[,1] == "Prostate Cancer"),1] = "Prostate cancer"
# allResults[which(allResults[,1] == "Major depressive disorder"),1] = "Major depression"
# allResults[which(allResults[,1] == "Lymphocytic leukemia"),1] = "Lymphoid leukemia"
# allResults[which(allResults[,1] == "Lymphoid leukemia, chronic"),1] = "Lymphoid leukemia"
# allResults[which(allResults[,1] == "Breast Cancer"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "Breast cancer (female)"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "ER-positive Breast Cancer"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "ER-negative Breast Cancer"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "Breast cancer intrinsic-like subtype (luminal A-like)"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "Breast cancer intrinsic-like subtype (luminal B/HER2-negative-like)"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "Breast cancer intrinsic-like subtype (luminal B-like)"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "Breast cancer intrinsic-like subtype (HER2-enriched-like)"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "Breast cancer intrinsic-like subtype (triple negative)"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "Estrogen receptor negative breast cancer"),1] = "Breast cancer"
# allResults[which(allResults[,1] == "estimated glomerular filtration rate (Cystatin C)"),1] = "estimated glomerular filtration rate"
# allResults[which(allResults[,1] == "Fasting glucose (body mass index adjusted)"),1] = "Fasting glucose"
# allResults[which(allResults[,1] == "Alzheimer's disease (late onset)"),1] = "Alzheimer's Disease"
# allResults[which(allResults[,1] == "Late-onset Alzheimer’s disease "),1] = "Alzheimer's Disease"
# allResults[which(allResults[,1] == "Colorectal cancer risk"),1] = "Colorectal cancer"
# allResults[which(allResults[,1] == "Rheumatoid Arthritis (CCP-positive)"),1] = "Rheumatoid Arthritis"
# allResults[which(allResults[,1] == "Rheumatoid Arthritis (CCP-negative)"),1] = "Rheumatoid Arthritis"
# allResults[which(allResults[,1] == "Basal cell carcinoma (MTAG)"),1] = "Basal cell carcinoma"
# allResults[which(allResults[,1] == "Coronary Artery disease"),1] = "Coronary artery disease"
# allResults[which(allResults[,1] == "estimated glomerular filtration rate (Creatine)"),1] = "estimated glomerular filtration rate"
# allResults[which(allResults[,1] == "estimated glomerular filtration rate (Cystatin C)"),1] = "estimated glomerular filtration rate"
# # actually don't do this, use the EFO IDs, as that is more precise!
# 




# eliminate PGS that are not LD-aware, (IE not LDpred or PRCS-CS)
ldAwareIndices = which(allResults[,5] == "DBSLMM" | allResults[,5] == "Integrative PGS (also referred to as metaGRS) of 22 component PGSs" | allResults[,5] == "lassosum" | allResults[,5] == "LDpred" | allResults[,5] == "LDpred2" | allResults[,5] == "PRS-CS" | allResults[,5] == "PRS-CS-auto" | allResults[,5] == "PRS-CSx")

LD_results = allResults[ldAwareIndices,]
nrow(LD_results) # 202



# exclude those that have too many or too few SNPs 
goodSampleSize = which(as.numeric(LD_results[,6]) >= 100 & as.numeric(LD_results[,6]) <= 1300000)
length(goodSampleSize) # 181

LD_N_results = LD_results[goodSampleSize,]


# filter for the single largest samplesize for each phenotype, to only keep one
tabeRes = table(LD_N_results[,8])
tableNames = names(tabeRes)
largestNRes = NULL
for (i in 1:length(tableNames)) {
  pheno = tableNames[i]
  largestNIndices = which(LD_N_results[,8] == pheno)
  phenoRes = LD_N_results[largestNIndices,]
  if ( length(largestNIndices) > 1) {
  largestSampleSizeIndices = which( as.numeric(phenoRes[,2]) == max( as.numeric(phenoRes[,2]) ) )
  phenoRes = phenoRes[largestSampleSizeIndices,]
  
  # in case we have the same sample size for multiple one, we make selection based on number of variants
  # we prefer the fewer ones
  if( length(largestSampleSizeIndices) > 1) {
    fewestSNPs =  which( as.numeric(phenoRes[,6]) == min( as.numeric(phenoRes[,6]) ) )
    phenoRes = phenoRes[fewestSNPs,]
    } 
  }
  largestNRes = rbind(largestNRes,phenoRes )
}

# there are still duplicates, filter based on 
duplicatedIdx = which( duplicated(largestNRes[,4]) )
#sub = largestNRes[duplicatedIdx,]

uniqueData = largestNRes[-c(duplicatedIdx), ]
nrow(uniqueData) # 48


write.table(uniqueData, paste0(outLoc, "_best.txt"), sep = "\t", row.names = F, col.names = T, quote = FALSE)


subset = uniqueData[,c(1,8)]

write.table(subset, paste0(outLoc, "_subset.csv"), sep = "\t", row.names = F, col.names = T, quote = FALSE)





