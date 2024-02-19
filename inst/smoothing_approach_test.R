library(data.table)
library(dplyr)
library(CINSignatureQuantification)
library(CINfits)

smoothSegmentsGEL <- function(segment.table,SMOOTHING_FACTOR){
    segment.table.smoothed.rounded <- segment.table %>%
        #mutate(segVal = round(segVal,digits = 2)) %>%
        group_by(chromosome,sample) %>%
        mutate(seg_diff = abs(segVal - lag(segVal))) %>%
        mutate(chng = ifelse(seg_diff > SMOOTHING_FACTOR,"TRUE","FALSE")) %>%
        mutate(chng = as.logical(ifelse(is.na(chng),"TRUE",chng))) %>%
        mutate(comb = cumsum(chng)) %>%
        group_by(chromosome,sample,comb) %>%
        select(-chng) %>%
        mutate(length = end - start) %>%
        filter(length > 60000) %>%
        summarise(across(start,min),across(end,max),across(segVal,median)) %>%
        #summarise(across(start,min),across(end,max),across(segVal,~weighted.mean(.,w=length,na.rm=TRUE))) %>%
        select(chromosome,start,end,segVal,sample) %>%
        mutate(chromosome = factor(chromosome,levels=c(1:22,"X","Y"))) %>%
        arrange(sample,chromosome,start)
    return(as.data.frame(segment.table.smoothed.rounded))
}

as <- read.table("inst/ascat_seg_data.tsv",header = T,sep = "\t")
asP <- processAlleleSpecific(as)

asP <- read.table("inst/britroc_30kb_ds_absCopyNumber_segmentTable.tsv",header = T,sep = "\t")

for(i in unique(asP$sample)){
    asPSample <- asP[asP$sample == i,which(!colnames(asP) %in% c("nAraw","nBraw"))]

    smoothSegments <- CINSignatureQuantification:::smoothSegments
    idSmoothingTargets <- CINSignatureQuantification:::idSmoothingTargets

    smoothSegmentsCINpkg <- function(data,WIGGLE){
        dfAllSegs <- idSmoothingTargets(data,WIGGLE = WIGGLE,colNameSegVal = "segVal",
                                        colNameChr = "chromosome", IGNOREDELS = FALSE)
        lRaw = split(dfAllSegs, dfAllSegs$sample)
        lSmooth = smoothSegments(lRaw, CORES=1, WIGGLE = WIGGLE,
                                 colNameMerge = "segVal", colNameChr = "chromosome",
                                 colNameStart = "start", colNameEnd = "end",
                                 IGNOREDELS = FALSE, asDf = FALSE)
        dtSmooth = rbindlist(lSmooth)
        return(as.data.frame(dtSmooth))
    }

    factor <- 0.25

    cinpkg <- smoothSegmentsCINpkg(asPSample,WIGGLE = factor)
    gelun <- smoothSegmentsGEL(asPSample,SMOOTHING_FACTOR = factor)

    print(nrow(cinpkg) - nrow(gelun))
}

