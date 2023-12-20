#' rescaleFit
#'
#' Refit a segmented copy number profile to a new ploidy-purity combination
#'
#' @param data data.frame containing a segmented copy number profile
#' @param ploidy new ploidy value to refit copy number profile
#' @param purity new ploidy value to refit copy number profile
#' @param segValOnly return only segVal vector rather than entire refitted segment
#'   table (Default: FALSE)
#'
#' @return data.frame
#' @export
#'
rescaleFit <- function(data=NULL,ploidy=NULL,purity=NULL,segValOnly=FALSE){

    relploidy <- calculatePloidy(data)
    cellploidy <- ploidy*purity + 2*(1-purity)
    scaledploidy <- relploidy/cellploidy

    scaled_cn_segVal <- depthtocn(data$segVal,purity,scaledploidy)
    data$segVal <- scaled_cn_segVal
    if(segValOnly){
        return(data$segVal)
    } else {
        return(data)
    }
}

## depthtocn
# support function to convert copy number given purity and single copy depth
depthtocn<-function(v,purity,scaledploidy)
{
    (v/scaledploidy-2*(1-purity))/purity
}

