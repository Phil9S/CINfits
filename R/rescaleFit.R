#' rescaleFit
#'
#' Refit a segmented copy number profile to a new ploidy-purity combination
#'
#' @param data data.frame containing a segmented copy number profile
#' @param old_ploidy old ploidy value from existing fit. Can be calculated from profile.
#' @param old_purity old purity value from existing fit
#' @param new_ploidy new ploidy value to refit copy number profile
#' @param new_purity new purity value to refit copy number profile
#' @param alleleSpecific rescale allele-specific copy number data (default: FALSE)
#'
#' @return data.frame
#' @export
#'
rescaleFit <- function(data=NULL,new_ploidy=NULL,new_purity=NULL,old_ploidy=NULL,old_purity=NULL,alleleSpecific=FALSE){

    if(is.null(data)){
        stop("no data")
    }
    if(is.null(old_purity)){
        stop("old purity must be provided")
    }
    if(is.null(new_purity)){
        stop("new purity must be provided")
    }
    if(is.null(new_ploidy)){
        stop("new ploidy must be provided")
    }
    if(is.null(old_ploidy)){
        old_ploidy <- calculatePloidy(data = data)
    }

    # cellploidy <- new_ploidy*new_purity + 2*(1-new_purity)
    # scaledploidy <- relploidy/cellploidy

    if(alleleSpecific){
        alleleRatioA <- data$nAraw / data$segVal
        alleleRatioA[is.na(alleleRatioA)] <- 0.5
        alleleRatioB <- 1-alleleRatioA
    }

    scaled_cn_segVal <- convertCN(data,
                                  old_ploidy,
                                  old_purity,
                                  new_ploidy,
                                  new_purity)

    data <- scaled_cn_segVal

    if(alleleSpecific){
        data$nAraw <- data$segVal * alleleRatioA
        data$nBraw <- data$segVal * alleleRatioB
    }

    return(data)
}

## depthtocn
# support function to convert copy number to new ploidy purity combination
convertCN<-function(v,old_ploidy,old_purity,new_ploidy,new_purity){

    dip <- 2
    original_segVal <- v$segVal
    old_sample_ploidy <- dip*(1-old_purity)+old_ploidy*old_purity
    seg_rel <- (original_segVal*old_purity)+(dip*(1-old_purity))

    c <- v
    c$segVal <- seg_rel
    new_rel_ploidy <- calculatePloidy(c)

    new_sample_ploidy <- dip*(1-new_purity)+new_ploidy*new_purity
    new_segdepth <- round(new_rel_ploidy/new_sample_ploidy,digits = 3)

    c$segVal <- (c$segVal/new_segdepth-dip*(1-new_purity))/new_purity

    #c$segVal[c$segVal < 0] <- 0
    #print(data.frame(v$segVal,c$segVal))
    return(c)
}

