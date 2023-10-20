#' calculateHomozygousLoss
#'
#' @param data
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
calculateHomozygousLoss <- function(data=NULL,threshold=0.2){
            data$length <- data$end - data$start
            hmzy <- sum(data$length[data$segVal < threshold]) / sum(data$length)*100
            return(hmzy)
}
