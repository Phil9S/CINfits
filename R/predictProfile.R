#' predictProfile
#'
#' @param qctable
#' @param model
#' @param flagThreshold
#'
#' @returns data.frame
#' @export
#'
predictProfile <- function(qctable = NULL,model = NULL,flagThreshold=0.74,multiFitFilter="clonality"){

    if(is.null(model)){
        selectedModel <- workflowsets::extract_workflow(final_rf_res)
    } else {
        selectedModel <- workflowsets::extract_workflow(model)
    }

    #print(selectedModel)
    #print(qcTablePred)
    newClass <- stats::predict(selectedModel,qctable,type = "class")
    newProb <- round(stats::predict(selectedModel,qctable,type = "prob"),digits = 3)

    qctable <- cbind(qctable,newClass,newProb)

    qctable <- triageProfile(qcTable = qctable,flagThreshold = flagThreshold,
                             multiFitFilter = multiFitFilter)

    return(qctable)
}

triageProfile <- function(qcTable = NULL,flagThreshold=0.74,multiFitFilter="clonality"){
    if(is.null(qcTable)){
        stop("no data")
    }

    qcTable <- qcTable %>%
            group_by(sample) %>%
            mutate(triageValue = abs(.pred_TRUE - .pred_FALSE)) %>%
            mutate(sumt = sum(as.logical(.pred_class)))

    switch(multiFitFilter,
           "clonality"={
               qcTable <- qcTable %>%
                   arrange(sample,clonality) %>%
                   mutate(use = ifelse(sumt <= 1,as.logical(.pred_class),
                                           ifelse(clonality == min(clonality),TRUE,FALSE)))
           },
           "segvariance"={
               qcTable <- qcTable %>%
                   arrange(sample,segvariance) %>%
                   mutate(use = ifelse(sumt <= 1,as.logical(.pred_class),
                                           ifelse(segvariance == min(segvariance),TRUE,FALSE)))
    })

    qcTable <- qcTable %>%
                mutate(flag = ifelse(triageValue < flagThreshold,"LOWPROBFLAG",NA)) %>%
                rename(pred_class = .pred_class,pred_TRUE = .pred_TRUE,pred_FALSE = .pred_FALSE) %>%
                mutate(use = pred_class) %>%
                mutate(notes = ifelse(is.na(flag),NA,flag))
    return(qcTable)
}
