## modify cell line ascat.sc and QC input for use

library(dplyr)
library(tidyr)

# Read cell line ASCAT.sc total CN data
cellSegs <- read.table("inst/cellLine_ascat_sc_fixed_purity_tCN_allSamples.tsv",
                       header = T,sep = "\t")

# Adjust column names
colnames(cellSegs) <- c("chromosome","start","end","rounded",
                        "segVal","logr","logr.sd","sample")

# Select columns
cellSegs <- cellSegs[,c(1:3,5,8)]

# Apply smoothing to adjacent segements with proximal segVal
cellSegsSmooth <- smoothProfile(cellSegs,smoothingFactor = 0.12)

# read QC and Fit data
cellQC <- read.table("inst/cellLine_fit_qc_table.tsv",
                     header = T,sep = "\t")
# cellFit <- read.table("inst/cellLine_ascat_sc_fixed_purity_profile_statistics.tsv"
#                       ,header = F,sep = "\t",col.names = c("sample","purity","ploidy"))
#
#
# cellQC <- cellQC %>%
#             left_join(.,cellFit,by = "sample") %>%
#             select(-notes,-plot)
# run fitInteractive() to generate compatible QC file and update USE column
