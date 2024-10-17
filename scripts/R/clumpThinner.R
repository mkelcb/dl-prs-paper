# Removes that all SNPs are at least 500 kb apart from each other of a Clumping result

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 
 
gwasLoc = args[1]
outputLoc= args[2]
minDistance = as.numeric(args[3])


# Debug vars
# gwasLoc =  "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/results/sims/GWAS_1_v0_Clump.clumps"
# outputLoc =   "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/PRS_GXE/scratch/dummy"
# minDistance= 500

# 1. load phenos
gwasRes= read.table(gwasLoc, header = F) 
gwasRes =na.omit(gwasRes)
minDistance=  minDistance * 1000 # supplied in KB, but we need it in bp

# constant, which are based off on the PLINK2 CLUMP file
#CHROM  POS     ID      P       TOTAL   NONSIG  S0.05   S0.01   S0.001  S0.0001 SP2
chromCol = 1
bpCol = 2
SNPcol= 3


# dont assume variants are ordered within each chromosome... we order them
# loop 1:22 chromosomes
gwasRes_sorted = NULL
i=1
for (i in 1:22) {
print(paste0("sorting chrom ", i))
# grab all SNPs for a chromosome
chromIndices = which(gwasRes[,chromCol] == i)
gwasRes_chrom = gwasRes[chromIndices,]

# sort by position, ascending within each chromosome
gwasRes_chrom <- gwasRes_chrom[order(gwasRes_chrom[,bpCol]),]

## add to ongoing one
gwasRes_sorted = rbind(gwasRes_sorted, gwasRes_chrom)
}

gwasRes = gwasRes_sorted

  

# loop variants from the 2nd
lastVariant_bp = gwasRes[1,bpCol]
lastVariant_chrom = gwasRes[1,chromCol]
i=2

variantsToExclude = c()
for (i in 2:10) {
#for (i in 2:nrow(gwasRes)) {
  # check if its same chromosome
  if (gwasRes[i,chromCol] == lastVariant_chrom) {
    dist = abs(lastVariant_bp - gwasRes[i,bpCol])
    # otherwise check if its distance to previous one is < minDistance
    if ( dist < minDistance) {
      # if yes add it to exclude list, and do NOT move the 'lastVariant' forward
      variantsToExclude = c(variantsToExclude, i)
      print(paste0("excluding variant ", i, " (chrom: ",gwasRes[i,chromCol], " / ", gwasRes[i,SNPcol],"), as distance was: ", dist, " (lastVariant_chrom: ",lastVariant_chrom, " / lastVariant_bp: ", lastVariant_bp, " / gwasRes[i,bpCol]: ", gwasRes[i,bpCol], ")"))
    } else {
      # if this next variant was far enough, then we nominate it to be the next 'lastVariant'
      lastVariant_bp = gwasRes[i,bpCol]
      lastVariant_chrom = gwasRes[i,chromCol]
    }
  } else {  # if it wasn't even the same chromosome then this variant will be the next lastVariant for sure
    lastVariant_bp = gwasRes[i,bpCol]
    lastVariant_chrom = gwasRes[i,chromCol]
  }
}

# V1       V2         V3            V4 V5 V6 V7 V8 V9 V10
# 1 10490881   rs481571 2.37850e-10 56  6 10 15  8  17
# 1 15344936   rs938252 5.42811e-09  6  1  0  0  1   4
# 1 18040004  rs2185324 2.14474e-05  4  1  1  1  0   1
# 1 18123332  rs4920385 1.08720e-27 14  0  0  0  2  12
# 1 18176173 rs11203219 4.02267e-05  1  0  1  0  0   0
# 1 20901688  rs4655223 3.01601e-05 10  5  2  1  1   1



# write all variants kept
if(length(variantsToExclude)> 0) {
  gwasRes_kept = gwasRes[-variantsToExclude,SNPcol]
  # write all variants excluded
  gwasRes_excluded = gwasRes[variantsToExclude,SNPcol]
  #write.table(gwasRes_excluded,file=paste0(outputLoc, "_excluded"), col.names = F, quote = F, row.names = F, sep = "\t")
  
} else {
  gwasRes_kept = gwasRes[,SNPcol]
}

#write.table(gwasRes,file=paste0(outputLoc,"_sorted"), col.names = F, quote = F, row.names = F, sep = "\t")

write.table(gwasRes_kept,file=outputLoc, col.names = F, quote = F, row.names = F, sep = "\t")

print(paste0("out of ",nrow(gwasRes), ", ", length(variantsToExclude), " were excluded, kept variants written to", outputLoc ))

