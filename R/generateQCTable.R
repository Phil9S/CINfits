generateQCTable <- function(data=NULL,metadata=NULL){
    if(is.null(data)){
        stop("no data")
    }
    if(is.null(metadata)){
        stop("no metadata")
    }

    if(is.null(data$sample)){
        stop("no sample field in segtable")
    }

    meta.cols <- c("sample","purity")

    if(!all(meta.cols %in% colnames(metadata))){
        stop("missing meta columns")
    }

    metadata.filt <- metadata[,meta.cols]
    stats <- rownamesToCol(calculateCINStats(data))

    qc.table <- merge(stats,metadata.filt,by="sample")
    qc.table$use <- rep(NA,times=nrow(qc.table))
    qc.table$notes <- rep("",times=nrow(qc.table))

    qc.table$clonality <- round(qc.table$clonality,digits = 3)
    qc.table$homozygousLoss <- round(qc.table$homozygousLoss,digits = 3)

    qc.table <- qc.table[,c("sample","segments","clonality","ploidy","purity",
                            "homozygousLoss","use","notes")]

    return(qc.table)
}

# Retains rownames from matrix outputs for export to tabular format
rownamesToCol <- function(x) {
    x <- as.data.frame(x)
    rn <- rownames(x)
    rownames(x) <- NULL
    y <- cbind(rn,x)
    colnames(y) <- c("sample",colnames(x))
    return(y)
}
