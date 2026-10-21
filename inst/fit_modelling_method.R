#load libraries
library(QDNAseqmod)
library(Biobase)
library(ggplot2)
library(stringr)


#getsegtable
getSegTable <- function(x){
    if(inherits(x,what = "QDNAseqCopyNumbers",which = F)){
        sn<-Biobase::assayDataElement(x,"segmented")
        fd <- Biobase::fData(x)
        use <- fd$use
        fdfiltfull<-fd[use,]
        sn<-sn[use,]
        sn <- as.data.frame(sn)
        if(ncol(sn) == 1){
            colnames(sn) <- colnames(x)
        }

        segTable<-c()
        for(s in colnames(sn)){
            for(c in unique(fdfiltfull$chromosome)){
                snfilt <- sn[fdfiltfull$chromosome == c,colnames(sn) == s]
                fdfilt <- fdfiltfull[fdfiltfull$chromosome == c,]
                sn.rle <- rle(snfilt)
                starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
                ends <- cumsum(sn.rle$lengths)

                segtmp <- do.call(rbind,lapply(1:length(sn.rle$lengths), function(s){
                    from <- fdfilt$start[starts[s]]
                    to <- fdfilt$end[ends[s]]
                    segValue <- sn.rle$value[s]
                    c(fdfilt$chromosome[starts[s]], from, to, segValue)
                }))

                segTableRaw <- cbind(segtmp,sample = rep(s,times=nrow(segtmp)))
                segTable <- rbind(segTable,segTableRaw)
            }
        }
        segTable <- as.data.frame(segTable)
        colnames(segTable) <- c("chromosome", "start", "end", "segVal","sample")
        segTable$segVal <- as.numeric(segTable$segVal)
        segTable$segVal <- round(as.numeric(segTable$segVal),3)
        segTable$segVal[segTable$segVal < 0] <- 0
        segTable$start <- as.numeric(segTable$start)
        segTable$end <- as.numeric(segTable$end)
        return(segTable)
    } else {
        # NON QDNASEQ BINNED DATA FUNCTION

    }
}

# converts read depth to copy number given purity and single copy depth
depthtocn <- function(x,purity,seqdepth){
    (x/seqdepth-2*(1-purity))/purity
}

calculateSegmentVar <- function(abs_seg = NULL, abs_cn = NULL, chromosomes = NULL) {
    abs_seg <- as.numeric(abs_seg)
    abs_cn <- as.numeric(abs_cn)

    if (is.null(chromosomes)) {
        chromosomes <- rep(1, length(abs_seg))
    }
    stopifnot(length(chromosomes) == length(abs_seg),
              length(abs_cn) == length(abs_seg))

    # a new segment starts where segment value or chromosome changes
    n <- length(abs_seg)
    newSeg <- c(TRUE, abs_seg[-1] != abs_seg[-n] | chromosomes[-1] != chromosomes[-n])
    newSeg[is.na(newSeg)] <- TRUE
    segId <- cumsum(newSeg)

    segVar <- tapply(abs_cn, segId, stats::var)

    medianVar <- stats::median(segVar, na.rm = TRUE)
    return(medianVar)
}

qc.data <- read.table("../../../OneDrive - CRUK Cambridge Institute/britroc-cn-analysis/absolute_PRE_down_sampling/fit_QC_predownsample.tsv",header = T,sep = "\t")
rdsdata <- paste0("../../../OneDrive - CRUK Cambridge Institute/britroc-cn-analysis/absolute_PRE_down_sampling/britroc_smoothed_copyNumbersSegmented.rds")

rds.rel <- readRDS(rdsdata)
samples <- qc.data[which(qc.data$SAMPLE_ID %in% colnames(rds.rel)),]

abs_profiles <- rds.rel[fData(rds.rel)$use,]

#rds.pdata <- pData(abs_profiles)
rds.obj <- abs_profiles[,colnames(abs_profiles) %in% samples$SAMPLE_ID]

fitstats <- do.call(rbind,lapply(colnames(rds.obj),FUN = function(sample){

    ploidies <- samples$ploidy[samples$SAMPLE_ID == sample]
    purities <- samples$purity[samples$SAMPLE_ID == sample]

    ind<-which(colnames(rds.obj)==sample)
    relcn<-rds.obj[,ind]
    # added by PS
    to_use <- fData(relcn)$use #
    relcn <- relcn[to_use,] #
    # added by PS
    copynumber<-assayDataElement(relcn,"copynumber")
    rel_ploidy<-mean(copynumber,na.rm=T)

    cn<-assayDataElement(relcn,"copynumber")
    seg<-assayDataElement(relcn,"segmented")
    cn<-as.numeric(cn[!is.na(cn),])
    seg<-as.numeric(seg[!is.na(seg),])

    res <- do.call(rbind,lapply(1:length(ploidies),FUN = function(i){
        ploidy <- ploidies[i]
        purity <- purities[i]

        cellploidy<-ploidy*purity+2*(1-purity)
        seqdepth<-rel_ploidy/cellploidy

        abs_seg <- depthtocn(seg,purity,seqdepth)
        abs_cnbin <- depthtocn(cn,purity,seqdepth)

        segmentVar <- calculateSegmentVar(abs_seg,abs_cnbin)

        relcn_temp <- relcn
        assayDataElement(relcn_temp,"segmented")[,1] <- abs_seg

        tab <- getSegTable(relcn_temp)

        stats <- as.data.frame(cbind(calculateCINStats(tab),segvariance=segmentVar))

        stats <- stats %>% tibble::rownames_to_column("sample") %>%
            dplyr::mutate(grid_ploidy = round(rep(ploidy,times=nrow(.)),digits = 1)) %>%
            dplyr::mutate(grid_purity = round(rep(purity,times=nrow(.)),digits = 2))

        return(stats)
    }))
}))

fitstats_ord <- fitstats[order(fitstats$sample,fitstats$clonality),]
all(fitstats_ord$sample == qc.data$SAMPLE_ID)
fitstats_ord$use <- qc.data$use

cellQC <- fitstats_ord %>% ## generated in
    tidyr::drop_na() %>%
    dplyr::select(-grid_ploidy) %>%
    dplyr::rename(purity = grid_purity) %>%
    dplyr::mutate(use = factor(use))

write.table(cellQC,file = "inst/britroc_multifit_meta.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
