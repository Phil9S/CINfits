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
        fd$use -> use
        fdfiltfull<-fd[use,]
        sn<-sn[use,]
        sn <- as.data.frame(sn)
        if(ncol(sn) == 1){
            colnames(sn) <- colnames(x)
        }

        segTable<-c()
        for(s in colnames(sn)){
            for(c in unique(fdfiltfull$chromosome))
            {
                snfilt <- sn[fdfiltfull$chromosome==c,colnames(sn) == s]
                fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
                sn.rle<-rle(snfilt)
                starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
                ends <- cumsum(sn.rle$lengths)
                lapply(1:length(sn.rle$lengths), function(s) {
                    from <- fdfilt$start[starts[s]]
                    to <- fdfilt$end[ends[s]]
                    segValue <- sn.rle$value[s]
                    c(fdfilt$chromosome[starts[s]], from, to, segValue)
                }) -> segtmp
                segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),
                                          sample = rep(s,times=nrow(matrix(unlist(segtmp), ncol=4, byrow=T))),stringsAsFactors=F)
                segTable<-rbind(segTable,segTableRaw)
            }
        }
        colnames(segTable) <- c("chromosome", "start", "end", "segVal","sample")
        segTable$segVal <- as.numeric(segTable$segVal)
        segTable$start <- as.numeric(segTable$start)
        segTable$end <- as.numeric(segTable$end)
        return(segTable)
    } else {
        # NON QDNASEQ BINNED DATA FUNCTION

    }
}

# convert depth to abs cn
depthtocn <- function(x,purity,seqdepth) #converts read depth to copy number given purity and single copy depth
{
    (x/seqdepth-2*(1-purity))/purity
}

calculateSegmentVar <- function(x,y){
    segRLE <- rle(x)

    segVar <- c()
    for(i in 1:length(segRLE$lengths)){
        if(i == 1){
            strt_idx <- 1
            end_idx <- segRLE$lengths[i]
            segVar <- append(segVar,var(y[strt_idx:end_idx]))
        } else {
            start_idx <- max(cumsum(segRLE$lengths[1:i-1]))
            strt_idx <- start_idx + 1
            if(i == length(segRLE$lengths)){
                end_idx <- segRLE$lengths[i] + strt_idx - 1
            } else {
                end_idx <- segRLE$lengths[i] + strt_idx
            }
            segVar <- append(segVar,var(y[strt_idx:end_idx]))
        }
    }

    medianVar <- median(segVar)
    return(medianVar)
}

qc.data <- read.table("../../OneDrive - CRUK Cambridge Institute/britroc-cn-analysis/absolute_PRE_down_sampling/fit_QC_predownsample.tsv",header = T,sep = "\t")
rdsdata <- paste0("../../OneDrive - CRUK Cambridge Institute/britroc-cn-analysis/absolute_PRE_down_sampling/britroc_smoothed_copyNumbersSegmented.rds")

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

        stats <- stats %>% rownames_to_column("sample") %>%
            mutate(grid_ploidy = round(rep(ploidy,times=nrow(.)),digits = 1)) %>%
            mutate(grid_purity = round(rep(purity,times=nrow(.)),digits = 2))

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

fitstats_ord_single <- fitstats_ord %>% ## remove sample duplicates and fit single profile
                        group_by(sample) %>%
                        filter(row_number() == 1 | use == TRUE) %>%
                        filter(row_number() == 1)
