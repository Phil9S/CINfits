#' calculateSegmentCounts
#'
#' Compute the number of segments for a given copy number profile
#' @param data copy number profile from which to count number of segments
#'
#' @return numeric
#' @export
#'
calculateSegmentCounts <- function(data=NULL){
    segCounts <- nrow(data)
    return(segCounts)
}
