#' plotSunrise
#'
#' Plot sunrise plot of clonality values computed across a grid search of
#' supplied ploidy and purity values
#'
#' @param data data.frame containing a segmented copy number profile
#' @param ploidys numeric vector of ploidy values (Default:
#'   "seq.int(1.2,8,0.1)")
#' @param purities numeric vector of purity values between 0 and 1.0 (Default:
#'   "seq.int(0.2,1,0.01)")
#' @param ploidy ploidy value of current profile from which to calculate sunrise
#' @param purity purity value of current profile from which to calculate sunrise
#' @param alleleSpecific use allele-specific data
#'
#' @return plot
#' @export
#'
plotSunrise <- function(data=NULL,ploidys=NULL,purities=NULL,ploidy=NULL,purity=NULL,alleleSpecifc=FALSE){
    if(is.null(ploidys)){
        ploidys <- seq.int(1.2,8,0.1)
    }
    if(is.null(purities)){
        purities <- seq.int(0.2,1,0.01)
    }
    clonality <- calculateSunrise(data=data,
                                  ploidys = ploidys,
                                  purities = purities,
                                  ploidy=ploidy,
                                  purity=ploidy,
                                  alleleSpecifc=alleleSpecifc)

    ## avoid R CMD check notes
    plpu <- data.frame(ploidy=c(as.character(ploidy)),
                       purity=c(as.character(purity)),
                       clonality=0.1)
    #print(plpu)
    #ploidy <- purity <- NULL

    ggplot2::ggplot(data = clonality,ggplot2::aes(x=ploidy,y=purity,fill=clonality)) +
        ggplot2::geom_tile() +
        #geom_point(data=plpu,shape=4, size = 10) +
        ggplot2::scale_x_discrete(breaks = unique(round(ploidys,digits = 0))) +
        ggplot2::scale_y_discrete(breaks = unique(round(purities,digits = 1))) +
        ggplot2::scale_fill_gradient(low = "blue",high = "yellow") +
        #ggplot2::scale_fill_gradient2(low = "blue",mid = "white",high = "yellow",midpoint = midpoint) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")
}

calculateSunrise <- function(data=NULL,ploidys=NULL,purities=NULL,ploidy=NULL,purity=NULL,alleleSpecifc=FALSE){
    if(is.null(ploidys)){
        ploidys <- seq.int(1.2,8,0.1)
    }
    if(is.null(purities)){
        purities <- seq.int(0.2,1,0.01)
    }
    comp_clonality <- function(x,y){
        calculateClonality(data = rescaleFit(data,
                                             old_ploidy = ploidy,
                                             old_purity = purity,
                                             new_ploidy = x,new_purity = y,
                                             alleleSpecific = alleleSpecifc))
    }

    clonality <- sapply(ploidys, function(x) mapply(FUN = comp_clonality,x,purities))

    rownames(clonality) <- purities
    colnames(clonality) <- ploidys

    clonality <- tibble::rownames_to_column(as.data.frame(clonality),var = "purity") %>%
        tidyr::pivot_longer(cols = -1,names_to = "ploidy",values_to = "clonality")
    return(clonality)
}
