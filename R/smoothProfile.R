#' smoothProfile
#'
#' @param data data.frame containing copy number segments
#' @param smoothingFactor Distance between two segments at which to collapse
#'   into a single segment (Default: 0.12)
#' @param implementation Smoothing implementation to use, either linear or
#'   iterative (Default: linear)
#' @param method method to collapse segments (Default: median). Options are
#'   median, mean, weighted.mean.
#' @param alleleSpecific use allele-specific implementation
#'
#' @return data.frame
#' @importFrom rlang .data
#' @export
#'
smoothProfile <- function(data=NULL,smoothingFactor=0.12,implementation="linear",method="median",alleleSpecific=FALSE){
    if(is.null(data)){
        stop("no data provided")
    }
    if(!is.numeric(smoothingFactor)){
        stop("smoothingFactor should be a float between 0 - 0.5")
    }
    if(!implementation %in% c("linear","iterative")){
        stop("smoothing implementation not available")
    }
    if(!method %in% c("median","mean","weighted.mean")){
        stop("smoothing method not available")
    }

    if(smoothingFactor > 0.5 | smoothingFactor < 0){
        stop("smoothingFactor should be between 0 and 0.5")
    }

    if(implementation == "iterative" & method != "weighted.mean"){
        warning("By default the iterative implementation uses weighted mean method for collapsing segments")
    }

    # drop 1bp segments
    data <- data[which(data$end - data$start > 1),]

    switch(implementation,
        "linear" = {
            dtSmooth <- linearSmooth(data,smoothingFactor,method,alleleSpecific)
        },
        "iterative"={
            dfAllSegs <- idSmoothingTargets(data,smoothingFactor = smoothingFactor)
            lRaw = split(dfAllSegs, dfAllSegs$sample)
            lSmooth = iterativeSmooth(lRaw,smoothingFactor = smoothingFactor)
            dtSmooth = data.table::rbindlist(lSmooth)
        })
    return(as.data.frame(dtSmooth))
}

linearSmooth <- function(data,smoothingFactor,method,alleleSpecific){
    if(alleleSpecific){
        segment.table <- data %>%
            dplyr::group_by(.data$chromosome,.data$sample) %>%
            dplyr::mutate(length = .data$end - .data$start) %>%
            dplyr::mutate(seg_diffA = abs(.data$nAraw - dplyr::lag(.data$nAraw))) %>%
            dplyr::mutate(seg_diffB = abs(.data$nBraw - dplyr::lag(.data$nBraw))) %>%
            ## If either alleles change greater than smoothingFactor mark as different
            dplyr::mutate(chng = ifelse(.data$seg_diffA > smoothingFactor | .data$seg_diffB > smoothingFactor,"TRUE","FALSE")) %>%
            dplyr::mutate(chng = as.logical(ifelse(is.na(.data$chng),"TRUE",.data$chng))) %>%
            dplyr::mutate(comb = cumsum(.data$chng)) %>%
            dplyr::group_by(.data$chromosome,.data$sample,.data$comb) %>%
            dplyr::select(-.data$chng)
        switch(method,
               "median"={
                   segment.table <- segment.table %>%
                       dplyr::summarise(dplyr::across(.data$start,min),dplyr::across(.data$end,max),
                                        dplyr::across(.data$segVal,stats::median),
                                        dplyr::across(.data$nAraw,stats::median),
                                        dplyr::across(.data$nBraw,stats::median))
               },
               "mean"={
                   segment.table <- segment.table %>%
                       dplyr::summarise(dplyr::across(.data$start,min),dplyr::across(.data$end,max),
                                        dplyr::across(.data$segVal,mean),
                                        dplyr::across(.data$nAraw,mean),
                                        dplyr::across(.data$nBraw,mean))
               },
               "weighted.mean"={
                   segment.table <- segment.table %>%
                       dplyr::summarise(dplyr::across(.data$start,min),dplyr::across(.data$end,max),
                                        dplyr::across(.data$segVal,~stats::weighted.mean(.,w=length,na.rm=TRUE)),
                                        dplyr::across(.data$nAraw,~stats::weighted.mean(.,w=length,na.rm=TRUE)),
                                        dplyr::across(.data$nBraw,~stats::weighted.mean(.,w=length,na.rm=TRUE)))
               })
    } else {
        segment.table <- data %>%
            dplyr::group_by(.data$chromosome,.data$sample) %>%
            dplyr::mutate(seg_diff = abs(.data$segVal - dplyr::lag(.data$segVal))) %>%
            dplyr::mutate(chng = ifelse(.data$seg_diff > smoothingFactor,"TRUE","FALSE")) %>%
            dplyr::mutate(chng = as.logical(ifelse(is.na(.data$chng),"TRUE",.data$chng))) %>%
            dplyr::mutate(comb = cumsum(.data$chng)) %>%
            dplyr::group_by(.data$chromosome,.data$sample,.data$comb) %>%
            dplyr::select(-.data$chng) %>%
            dplyr::mutate(length = .data$end - .data$start)
        switch(method,
           "median"={
               segment.table <- segment.table %>%
                   dplyr::summarise(dplyr::across(.data$start,min),dplyr::across(.data$end,max),
                                    dplyr::across(.data$segVal,stats::median))
           },
           "mean"={
               segment.table <- segment.table %>%
                   dplyr::summarise(dplyr::across(.data$start,min),dplyr::across(.data$end,max),
                                    dplyr::across(.data$segVal,mean))
           },
           "weighted"={
               segment.table <- segment.table %>%
                   dplyr::summarise(dplyr::across(.data$start,min),dplyr::across(.data$end,max),
                                    dplyr::across(.data$segVal,~stats::weighted.mean(.,w=length,na.rm=TRUE)))
           })
    }

    segment.table <- segment.table %>%
        dplyr::select(-.data$comb) %>%
        dplyr::relocate(.data$sample,.after = dplyr::last_col()) %>%
        dplyr::relocate(.data$segVal,.before = .data$sample) %>%
        dplyr::mutate(chromosome = factor(.data$chromosome,levels=c(1:22,"X","Y"))) %>%
        dplyr::arrange(.data$sample,.data$chromosome,.data$start)
    return(as.data.frame(segment.table))
}

idSmoothingTargets <- function(dfAllSegs,smoothingFactor) {
    ## smoothing procedure implemented by Ruben Drews in Nature 2022 modified by Phil Smith
    colNameSegVal <- "segVal"
    colNameChr <- "chromosome"
    # Take differences to segment down below
    dfAllSegs$diffs = c( abs( dfAllSegs[[colNameSegVal]][1:(nrow(dfAllSegs)-1)] - dfAllSegs[[colNameSegVal]][2:nrow(dfAllSegs)] ), smoothingFactor+1)
    # Set TRUE if difference to next segment is smaller than the user supplied cutoff
    dfAllSegs$smooth = dfAllSegs$diffs <= smoothingFactor
    # Set all segments which are last in a chromosome to FALSE. This also prevents leaking to other samples and cohorts.
    dfAllSegs$smooth[ cumsum( rle(as.character(dfAllSegs[[colNameChr]]))$lengths ) ] = FALSE
    return(dfAllSegs)
}

iterativeSmooth <- function(lRaw, smoothingFactor){
    # smoothing procedure implemented by Ruben Drews in Nature 2022 modified by Phil Smith
    # Add diff column to names we want to keep when merging (comes from function "idSmoothingTargets").
    colNameMerge = "segVal"
    colNameChr = "chromosome"
    colNameStart = "start"
    colNameEnd = "end"
    colNameMerge <-  c(colNameMerge, "diffs")
    lSmooth <- lapply(lRaw,FUN = function(x){
        thisOut = x
        stillSmoothing = sum(thisOut$smooth)
        while( stillSmoothing > 0 ) {
            thisSample = thisOut

            rleRaw = rle(thisSample$smooth)
            indRaw = cumsum(rleRaw$lengths)[ ! rleRaw$values ] + 1
            indRaw = indRaw[ -length(indRaw) ]
            if( rleRaw$values[1] ) { indRaw = c(1, indRaw) }

            # loop over start indices of TRUE chains.
            for(i in indRaw) {
                # detect length of segments to smooth. add 1 as the last segment has a FALSE value in it but still belongs to this chain.
                endOfStreak = i + rle(thisSample$smooth[i:nrow(thisSample)])$lengths[1]
                # extract reads
                dfMerge = thisSample[i:endOfStreak,]
                # too stupid to make this work with data.table
                newElement = as.data.frame( dfMerge[1,] )
                # Get new end and check first wether valid number.
                newEnd = dfMerge[nrow(dfMerge),][[colNameEnd]]
                if(! is.null(newEnd)) {
                    newElement[[colNameEnd]] = newEnd
                } else {
                    stop("New end coordinate is null. Supplied correct column name?")
                }
                ## Column "segVal" will be dealt with in a minute. Column "diffs" later when running again idSmoothingTargets.

                # Merge cn specifically by taking the length of the elements into consideration
                widthWeights = dfMerge[[colNameEnd]] - dfMerge[[colNameStart]]
                newElement[[colNameMerge[1]]] = stats::weighted.mean(dfMerge[[colNameMerge[1]]], widthWeights)
                # Replace all to merge segments with the new merged segment. Later delete duplicated.
                thisOut[i:endOfStreak,] = newElement
            }
        # as we have replaced all segments with the new mean segment, we need to remove the duplicates
        thisOut = thisOut[ ! duplicated(thisOut), ]
        # again detect segments which needs smoothing
        thisOut = idSmoothingTargets(thisOut, smoothingFactor)
        stillSmoothing = sum(thisOut$smooth)
        }
        # after smoothing is finished, change name of cohort
        thisOut$smooth = NULL
        thisOut$diffs = NULL
        return( thisOut )
    })
    return(lSmooth)
}
