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
