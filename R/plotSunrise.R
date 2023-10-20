plotSunrise <- function(data=NULL,ploidys=NULL,purities=NULL){
    if(is.null(ploidys)){
        ploidys <- seq.int(1.2,8,0.1)
    }
    if(is.null(purities)){
        purities <- seq.int(0.2,1,0.01)
    }

    comp_clonality <- function(x,y){
        calculateClonality(data = rescale_fit(data,ploidy = x,purity = y))
    }

    clonality <- sapply(ploidys, function(x) mapply(FUN = comp_clonality,x,purities))

    heatmap(clonality,
            Rowv = NA,Colv = NA,
            labRow = purities,labCol = ploidys,
            xlab = "ploidy",ylab = "purity")
}
