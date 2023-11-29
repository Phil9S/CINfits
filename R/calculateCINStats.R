## calculateCINStats

calculateCINStats <- function(data=NULL){
    list.tabs <- split(data,f=data$sample)

    clonality <- unlist(lapply(list.tabs,calculateClonality))
    ploidy <-  unlist(lapply(list.tabs,calculatePloidy))
    homozygousLoss <-  unlist(lapply(list.tabs,calculateHomozygousLoss))

    comb <- cbind(clonality,ploidy,homozygousLoss)
    return(comb)
}
