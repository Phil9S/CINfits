#library(Biobase)

depthtocn<-function(x,purity,seqdepth) #converts read depth to copy number given purity and single copy depth
{
    (x/seqdepth-2*(1-purity))/purity
}

cntodepth<-function(cn,purity,seqdepth) #converts copy number to read depth given purity and single copy depth
{
    seqdepth*((1-purity)*2+purity*cn)
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

