calculateInterleaving <- function(data=NULL,states=4,distance=0.2){
    if(is.null(data)){
        stop("no data")
    }

    data.filt <- data %>%
                    dplyr::mutate(close_distance =
                        ifelse(abs(.data$segVal - round(.data$segVal,digits = 0)) < distance,TRUE,FALSE)) %>%
                    dplyr::mutate(segValDist = ifelse(.data$close_distance,round(.data$segVal,digits = 0),NA)) %>%
                    dplyr::filter(.data$close_distance,.data$segValDist <= states,.data$segVal > distance) %>%
                    dplyr::mutate(length = .data$end - .data$start) %>%
                    #dplyr::group_by(.data$chromosome) %>%
                    dplyr::mutate(pctLen = .data$length / sum(.data$length)) %>%
                    dplyr::group_by(.data$segValDist) %>%
                    dplyr::summarise(across(.data$pctLen,sum))

    data.filt <- as.data.frame(data.filt)
    if(nrow(data.filt) < states){
        missingRows <- which(!c(1:states) %in% as.numeric(data.filt$segValDist))

        data.filt <- rbind(data.filt,
                           data.frame(segValDist=missingRows,
                                      pctLen=rep(0,times=length(missingRows))))
    }

    rownames(data.filt) <- paste0("cnstate_",data.filt$segValDist)
    data.format <- as.data.frame(t(as.matrix(data.filt)[,2]))
    #print(data.format)

    return(data.format)
}
