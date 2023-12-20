#' calculateCINStats
#'
#' calculate CIN statistics and metrics over a set of segmented copy number
#' profiles
#'
#' @param data data.frame or list of segmented copy number profiles
#'
#' @return data.frame
#' @export
#'
calculateCINStats <- function(data=NULL){
    if(is.null(data)){
        stop("no data")
    }
    if(is.data.frame(data)){
        list.tabs <- split(data,f=data$sample)
    } else {
        list.tabs <- data
    }

    clonality <- unlist(lapply(list.tabs,calculateClonality))
    ploidy <-  unlist(lapply(list.tabs,calculatePloidy))
    homozygousLoss <-  unlist(lapply(list.tabs,calculateHomozygousLoss))

    comb <- cbind(clonality,ploidy,homozygousLoss)
    return(comb)
}
