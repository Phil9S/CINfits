#' calculateClonality
#'
#' Compute clonality (measure of segment noise) using either weighted or
#' non-weighted mean of absolute segment distance from integer
#'
#' @param data data.frame containing a segmented copy number profile
#' @param weighted boolean to use weighted.mean to compute clonality (default:
#'   TRUE)
#'
#' @return numeric
#' @export
#'
calculateClonality <- function(data=NULL,weighted=TRUE){
    segvals <- data$segVal
    rounded <- round(segvals,digits = 0)
    if(is.null(data$length)){
        lengths <- data$end - data$start
    } else {
        lengths <- data$length
    }
    diff <- abs(rounded - segvals)
    if(weighted){
        clonality <- stats::weighted.mean(diff,lengths)
    } else {
        clonality <- mean(diff)
    }
    return(clonality)

}
