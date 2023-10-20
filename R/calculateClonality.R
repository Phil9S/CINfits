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
        clonality <- weighted.mean(diff,lengths)
    } else {
        clonality <- mean(diff)
    }
    return(clonality)

}
