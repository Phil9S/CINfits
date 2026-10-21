#' calculateHomozygousLoss
#'
#' Compute basepairs of genome which is homozygous loss across a segmented copy
#' number profile. The specified threshold is used to determine which segments
#' are considered homozygous loss.
#'
#' @param data data.frame containing a segmented copy number profile
#' @param threshold Copy number value at which to consider a segment as
#'   homozygous loss (Default: 0.4)
#'
#' @return numeric
#' @export
#'
calculateHomozygousLoss <- function(data=NULL,threshold=0.4){
            data$length <- data$end - data$start
            # to match gridsearch in swgs-abscn
            hmzy <- sum(data$length[data$segVal < threshold])
            return(hmzy)
}
