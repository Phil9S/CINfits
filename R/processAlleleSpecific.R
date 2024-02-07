#' processAlleleSpecific
#'
#' @param data
#' @param nameMaj
#' @param nameMinor
#' @param chrName
#' @param startName
#' @param endName
#' @param sampleName
#' @param returnTotal description
#' @param returnSmoothed
#' @param smoothingFactor
#'
#' @return data.frame
#' @export
#'
#' @examples

processAlleleSpecific <- function(data = NULL,
                                  nameMaj="nAraw",
                                  nameMinor="nBraw",
                                  chrName="chr",
                                  startName="startpos",
                                  endName="endpos",
                                  sampleName="sample",
                                  returnTotal=TRUE,
                                  returnSmoothed=FALSE,
                                  smoothingFactor=0.1){
    if(is.null(data)){
        stop("no data")
    }

    origCols <- c(sampleName,chrName,startName,endName,nameMaj,nameMinor)

    if(!all(origCols %in% colnames(data))){
        print(colnames(data))
        print(origCols)
        stop("mismatched col names")
    }

    dataSub <- data[,origCols]
    colName <- c("sample","chromosome","start","end","nAraw","nBraw")
    colnames(dataSub) <- colName

    colNameOrdered  <- c("chromosome","start","end","nAraw","nBraw")

    if(returnTotal){
        dataSub$segVal <- dataSub$nAraw + dataSub$nBraw
        dataSub <- dataSub[,c(colNameOrdered,"segVal","sample")]
    } else {
        dataSub <- dataSub[,c(colNameOrdered,"sample")]
    }

    if(returnSmoothed){
        dataSub <- dataSub[,which(!colnames(dataSub) %in% c("nAraw","nBraw"))]
        warning("smoothing not implemented")
    }

    return(dataSub)
}
