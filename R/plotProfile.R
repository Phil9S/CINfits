#' plotProfile
#'
#' Produce a plot of segments from a given segmented copy number profile
#'
#' @param data data.frame containing a segmented copy number profile
#' @param sample vector of length 1 containing either a sample name or sample index
#' @param cn.max maximum copy number to plot - Values greater than this are truncated to the specified value
#' @param purity given purity for profile being plotted - not required for non-interactive use.
#' @param alleleSpecific Boolean as to whether plot total or allele-specific copy number.
#' @param cols colours for segments plotted on copy number profile
#'
#' @return plot
#' @export
#'
plotProfile <- function(data=NULL,sample=NULL,cn.max=15,purity=NULL,alleleSpecific=FALSE,cols=NULL){

    segTab <- data
    segTab$chromosome <- factor(segTab$chromosome,
                                levels = stringr::str_sort(unique(segTab$chromosome),
                                                           numeric = T))
    ob.pl <- calculatePloidy(segTab)

    if(!is.null(purity)){
        ob.pu <- purity
    } else {
        ob.pu <- NA
    }

    if(alleleSpecific){
        if(is.null(data$segVal)){
            data$segVal <- data$nAraw + data$nBraw
        }
    }

    if(max(segTab$segVal) > cn.max){
        segTab$segVal[segTab$segVal > cn.max] <- cn.max
        ylim <- c(0,cn.max)
    } else {
        ylim <- c(0,round(max(segTab$segVal))+1)
    }

    if(is.null(cols)){
        if(alleleSpecific){
            nAcol <- "red"
            nBcol <- "blue"
        } else {
            col <- "blue"
        }
    } else {
        if(alleleSpecific){
            if(length(cols) == 2){
                nAcol <- cols[1]
                nBcol <- cols[2]
            } else {
                stop("incorrect len for cols")
            }
        } else {
            if(length(cols) == 1){
                col = cols
            } else {
                stop("incorrect len for cols")
            }
        }
    }

    seg.n <- nrow(segTab)

    chrom.len <- data.frame(Group.1=unique(segTab$chromosome))
    chrom.len$x.max <- stats::aggregate(segTab$end,
                                        by = list(segTab$chromosome),FUN = max)$x
    chrom.len$x.min <- stats::aggregate(segTab$start,
                                        by = list(segTab$chromosome),FUN = min)$x
    chrom.len$Group.1 <- factor(chrom.len$Group.1,
                                levels = stringr::str_sort(unique(chrom.len$Group.1),
                                                           numeric = T))
    chrom.len <- chrom.len[order(chrom.len$Group.1),]

    segTab$startf <- getcoordinates(chr = segTab$chromosome,
                                    pos = segTab$start,
                                    chrom.len = chrom.len)
    segTab$endf <- getcoordinates(chr = segTab$chromosome,
                                  pos = segTab$end,
                                  chrom.len = chrom.len)
    chrom.len$flatm <- getcoordinates(chr = chrom.len$Group.1,
                                      pos = chrom.len$x.min,
                                      chrom.len = chrom.len)
    chrom.len$flats <- getcoordinates(chr = chrom.len$Group.1,
                                      pos = chrom.len$x.max/2,
                                      chrom.len = chrom.len)
    chrom.len$flate <- getcoordinates(chr = chrom.len$Group.1,
                                      pos = chrom.len$x.max,
                                      chrom.len = chrom.len)

    title <- sample
    if(is.null(purity)){
        sub.title <- paste0("ploidy: ",ob.pl," | segments: ",seg.n)
    } else {
        sub.title <- paste0("ploidy: ",ob.pl," | purity: ",ob.pu," | segments: ",seg.n)
    }

    rect.col <- ifelse(seq_along(chrom.len$Group.1) %% 2 == 0,"white","grey95")

    graphics::par(mar=c(5, 4, 4, 4) + 0.2,xpd=FALSE)
    graphics::plot(NA,
                   xlab="chromosome",
                   ylab="absolute copy number",
                   las=1,
                   xlim=c(min(chrom.len$flatm),
                          max(chrom.len$flate)),
                   ylim=ylim,
                   xaxs="i",
                   xaxt="n",
                   yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]),
                   yaxs="i")
    graphics::rect(xleft = chrom.len$flatm,
                   xright = chrom.len$flate,
                   ybottom = 0,ytop = cn.max,
                   col=rect.col,
                   border = NA)
    graphics::axis(1, at=chrom.len$flats, labels=chrom.len$Group.1)
    graphics::box()
    graphics::mtext(side=3, line=2, at=-0.07, adj=0, cex=1.2, title)
    graphics::mtext(side=3, line=1, at=-0.07, adj=0, cex=1, sub.title)
    graphics::abline(h = seq.int(1,cn.max-1,1),lty="dashed",col="gray50")

    if(alleleSpecific){
        graphics::segments(x0 = segTab$startf,y0 = segTab$nAraw,
                           x1 = segTab$endf,y1 = segTab$nAraw,lwd=3,col=nAcol)
        graphics::segments(x0 = segTab$startf,y0 = segTab$nBraw,
                           x1 = segTab$endf,y1 = segTab$nBraw,lwd=3,col=nBcol)
    } else {
        graphics::segments(x0 = segTab$startf,y0 = segTab$segVal,
                           x1 = segTab$endf,y1 = segTab$segVal,lwd=3,col=col)
    }
}

## helper function to get coordinates from segment data
getcoordinates <- function(chr, pos, chrom.len) {
    posflat <- pos
    offset <- 0
    for (contig_ix in 1:nrow(chrom.len)) {
        on_contig <- chr == chrom.len$Group.1[contig_ix]
        posflat[on_contig] <- pos[on_contig] + offset
        offset <- offset + chrom.len$x.max[contig_ix]
    }
    posflat
}
