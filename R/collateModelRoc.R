collateModelRoc <- function(models,usePR=FALSE){
    if(is.null(models)){
        stop("no models provided")
    }

    if(!inherits(models,what = "list")){
        models <- list(models)
    }

    model_roc <- do.call(rbind,lapply(models,FUN = function(model){

        modelClass <- class(extract_spec_parsnip(model))[1]

        if(usePR){
            roccurve <- model %>%
                collect_predictions()  %>%
                pr_curve(use, .pred_FALSE) %>%
                mutate(model = modelClass)
        } else {
            roccurve <- model %>%
                collect_predictions()  %>%
                roc_curve(use, .pred_FALSE) %>%
                mutate(model = modelClass)
        }

    }))

    return(model_roc)
}
