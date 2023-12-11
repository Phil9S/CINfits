plotSunrise <- function(data=NULL,ploidys=NULL,purities=NULL){

    clonality <- calculateSunrise(data=data,ploidys = ploidys,purities = purities)

    ggplot2::ggplot(data = clonality,ggplot2::aes(ploidy,purity,fill=clonality)) +
        ggplot2::geom_tile() +
        ggplot2::theme_bw() +
        ggplot2::theme()
    # heatmap(clonality,
    #         Rowv = NA,Colv = NA,
    #         labRow = purities,labCol = ploidys,
    #         xlab = "ploidy",ylab = "purity")
}

calculateSunrise <- function(data=NULL,ploidys=NULL,purities=NULL){
    if(is.null(ploidys)){
        ploidys <- seq.int(1.2,8,0.1)
    }
    if(is.null(purities)){
        purities <- seq.int(0.2,1,0.01)
    }
    comp_clonality <- function(x,y){
        calculateClonality(data = rescaleFit(data,ploidy = x,purity = y))
    }

    clonality <- sapply(ploidys, function(x) mapply(FUN = comp_clonality,x,purities))

    rownames(clonality) <- purities
    colnames(clonality) <- ploidys

    clonality <- tibble::rownames_to_column(as.data.frame(clonality),var = "purity") %>%
        tidyr::pivot_longer(cols = -1,names_to = "ploidy",values_to = "clonality")
    return(clonality)
}
