collateModelRoc <- function(models,usePR=FALSE){
    use=`.pred_FALSE`=NULL
    if(is.null(models)){
        stop("no models provided")
    }

    if(!inherits(models,what = "list")){
        models <- list(models)
    }

    model_roc <- do.call(rbind,lapply(models,FUN = function(model){

        modelClass <- class(workflowsets::extract_spec_parsnip(model))[1]

        if(usePR){
            roccurve <- model %>%
                tune::collect_predictions()  %>%
                yardstick::pr_curve(use, .pred_FALSE) %>%
                dplyr::mutate(model = modelClass)
        } else {
            roccurve <- model %>%
                tune::collect_predictions()  %>%
                yardstick::roc_curve(use, .pred_FALSE) %>%
                dplyr::mutate(model = modelClass)
        }

    }))

    return(model_roc)
}
