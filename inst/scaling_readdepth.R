#library(Biobase)

depthtocn<-function(x,purity,seqdepth) #converts read depth to copy number given purity and single copy depth
{
    (x/seqdepth-2*(1-purity))/purity
}

cntodepth<-function(x,purity,seqdepth) #converts copy number to read depth given purity and single copy depth
{
    seqdepth*((1-purity)*2+purity*x)
}

#rds.obj <- readRDS("../../Downloads/sc_pseudo_30kb_relSmoothedCN.rds")
#relcn<-rds.obj[,1]
# added by PS
# to_use <- fData(relcn)$use #
# relcn <- relcn[to_use,] #
# added by PS
# copynumber<-assayDataElement(relcn,"copynumber")
# rel_ploidy<-mean(copynumber,na.rm=T)
# num_reads<-sum(copynumber,na.rm=T)
# sample <- sampleNames(relcn)

purity <- 0.7
ploidy <- 3

# Original equation
#downsample_depth <- (((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15*(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)

nbins_ref_genome <- 30000
read_ratio <- 1/0.91
rpc <- 15
scaling_factor <- rpc*nbins_ref_genome*read_ratio
sampleploidy <- 2*(1-purity)+ploidy*purity
seqdepth<-rel_ploidy/sampleploidy

downsample_depth <- (sampleploidy^2)/(ploidy*purity^2)*scaling_factor

# gene_bin_seg <- 15000
# cn<-assayDataElement(relcn,"copynumber")
# seg<-assayDataElement(relcn,"segmented")
# cn<-as.numeric(cn[!is.na(cn),])
# seg<-as.numeric(seg[!is.na(seg),])
# integer_cn<-round(depthtocn(seg,purity,seqdepth))
# abs_cn<-depthtocn(seg,purity,seqdepth)
# diffs<-abs(abs_cn-integer_cn)
# TP53cn<-round(depthtocn(median(seg[gene_bin_seg]),purity,seqdepth),1) # to 1 decimal place
# expected_TP53_AF<-TP53cn*purity/(TP53cn*purity+2*(1-purity))



####
br <- readRDS("../../OneDrive - CRUK Cambridge Institute/britroc-cn-analysis/absolute_PRE_down_sampling/britroc_smoothed_copyNumbersSegmented.rds")
t <- CINSignatureQuantification:::getSegTable(br[,1:20])
rd_profile <- t[t$sample == "IM_11",]
rd_profile[,2:4] <- apply(rd_profile[,2:4],MARGIN = 2,as.numeric)

data <- read.table("inst/britroc_30kb_ds_absCopyNumber_segmentTable.tsv",header = T,sep = "\t")
test_data <- data[data$sample == "IM_11",]

ploidy <- 2.97
purity <- 0.44

t1[1:10]
t2[2:11]
# existing profile
rel_ploidy <- calculatePloidy(rd_profile)
# new ploidy purity
sampleploidy <- 2*(1-purity)+ploidy*purity
#ratio
seqdepth<-rel_ploidy/sampleploidy
rel_to_abs <- (rd_profile$segVal/seqdepth-2*(1-purity))/purity

#######

data <- read.table("inst/britroc_30kb_ds_absCopyNumber_segmentTable.tsv",header = T,sep = "\t")
test_data <- data[data$sample == "IM_11",]

old_ploidy <- 2.97
old_purity <- 0.44

new_ploidy <- 1.4
new_purity <- 0.9

old_sample_ploidy <- 2*(1-old_purity)+old_ploidy*old_purity

original_segVal <- test_data$segVal
test_rel <- test_data
test_rel$segVal <- (original_segVal*0.44)+(2*(1-old_purity))

test_new <- test_data
new_rel_ploidy <- calculatePloidy(test_rel)
new_sample_ploidy <- 2*(1-new_purity)+new_ploidy*new_purity
new_segdepth <- round(new_rel_ploidy/new_sample_ploidy,digits = 3)
test_new$segVal <- (test_rel$segVal/new_segdepth-2*(1-new_purity))/new_purity

data.frame(test_data$segVal,test_new$segVal)
