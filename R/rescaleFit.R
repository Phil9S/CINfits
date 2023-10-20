rescaleFit <- function(data=NULL,ploidy=NULL,purity=NULL,segValOnly=FALSE){

    relploidy <- calculatePloidy(data)
    cellploidy <- ploidy*purity + 2*(1-purity)
    scaledploidy <- relploidy/cellploidy

    scaled_cn_segVal <- depthtocn(data$segVal,purity,scaledploidy)
    data$segVal <- scaled_cn_segVal
    if(segValOnly){
        return(data$segVal)
    } else {
        return(data)
    }
}

depthtocn<-function(v,purity,scaledploidy) #converts copy number given purity and single copy depth
{
    (v/scaledploidy-2*(1-purity))/purity
}

