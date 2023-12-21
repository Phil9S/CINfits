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
rescaleFit <- function(data=NULL,ploidy=NULL,new_purity=NULL,segValOnly=FALSE){

    old_purity <- 0.44 # TEST VALUE FOR IM_11

    relploidy <- calculatePloidy(data = data)
    cellploidy <- ploidy*new_purity + 2*(1-new_purity)
    scaledploidy <- relploidy/cellploidy

    scaled_cn_segVal <- convertCN(data$segVal,old_purity,new_purity,scaledploidy)
    data$segVal <- scaled_cn_segVal
    if(segValOnly){
        return(data$segVal)
    } else {
        return(data)
    }
}

## depthtocn
# support function to convert copy number given purity and single copy depth
convertCN<-function(v,old_purity,new_purity,scaledploidy)
{
    #v <- scaledploidy*((1-old_purity)*2+old_purity*v)
    c <- (v/scaledploidy-2*(1-new_purity))/new_purity
    return(c)
}

