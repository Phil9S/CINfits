#' calculatePloidy
#'
#' Compute ploidy using mean copy number segments, weighted by segment length.
#'
#' @param data data.frame containing a segmented copy number profile
#'
#' @return numeric
#' @export
#'
calculatePloidy <- function(data){
    segs <- data$segVal
    if(is.null(data$length)){
        lengths <- data$end - data$start
    } else {
        lengths <- data$length
    }
    p <- round(weighted.mean(segs,lengths),digits = 3)
    return(p)
}
