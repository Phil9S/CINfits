#' applyModel
#'
#' @param model fitted tidymodel workflow object
#' @param newData new data matching fitted data format
#'
#' @return prediction of new data using given model
#' @export
#'
applyModel <- function(model = NULL,newData = NULL){
    extractedWF <- workflowsets::extract_workflow(model)
    newPred <- predict(extractedWF,newData)
    return(newPred)
}
