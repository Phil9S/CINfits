#' calculateHomozygousLoss
#'
#' Compute proportion of genome which is homozygous loss across a segmented copy
#' number profile. The specified threshold is used to determine which segments
#' are considered homozygous loss.
#'
#' @param data data.frame containing a segmented copy number profile
#' @param threshold Copy number value at which to consider a segment as
#'   homozygous loss (Default: 0.75)
#'
#' @return numeric
#' @export
#'
calculateHomozygousLoss <- function(data=NULL,threshold=0.75){
            data$length <- data$end - data$start
            hmzy <- sum(data$length[data$segVal < threshold]) / sum(data$length)*100
            return(hmzy)
}
