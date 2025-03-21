
# Code repository for "Performance of deep-learning based approaches to improve polygenic scores"

DOI: 

This respository represents the last snapshot of the bash and R scripts used to generate our results and is provided as-is. As the analysis involved a lot of input/output operations on very large files, these were generated asynchronously on a cluster. Thus these scripts are meant to be executed on the command line manually, block-by-block, waiting for the remote jobs to finish and verifying the integrity of the resulting files at each step. To reduce code duplication, certain functions that were reused multiple times are defined only once across all files, however, they may be called from different scripts.

**The simulation analyses:**
1. PRS_GXE_sims.sh: code relevant for the simulation analyses

**The real data analyses:**
1. PRS_GXE_v2.sh: all real data analyses
2. PRS_GX_functions_v2.sh: common functions

**neural-network pytorch model:**
1. scripts/python/Knet/: Visual studio project that implements the neural-network models used in the study