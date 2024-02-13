#' processAlleleSpecific
#'
#' @param data absolute copy nubmer profile containing allele specific copy states.
#' @param nameMaj column name of major allele column (default: nAraw).
#' @param nameMinor column name of minor allele column (default: nBraw).
#' @param chrName column name of chromosome column (default: chr).
#' @param startName column name of genome start position column (default: startpos).
#' @param endName column name of genome end position column (default: endpos).
#' @param sampleName column name of sample identifier column (default: sample).
#' @param returnTotal A boolean value as to whether to additionally return total copy number, returned as 'segVal' (default: TRUE).
#' @param returnSmoothed A boolean value as to whether to return smoothed total copy number and drop allele-specific information.
#' @param smoothingFactor A float between 0 and 1 which is used to smooth and collapse adjacent total copy number segments (default: 0.1).
#'
#' @return data.frame
#' @export
#'

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
        if(!is.numeric(smoothingFactor)){
            stop("smoothingFactor not numeric")
        }
        if(smoothingFactor < 0 | smoothingFactor > 1){
            stop("smoothingFactor out of range")
        }
        dataSub <- dataSub[,which(!colnames(dataSub) %in% c("nAraw","nBraw"))]

        warning("smoothing not implemented")
    }

    return(dataSub)
}
