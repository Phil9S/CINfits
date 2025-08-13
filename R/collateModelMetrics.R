collateModelMetrics <- function(models){
    if(is.null(models)){
        stop("no models provided")
    }

    if(!inherits(models,what = "list")){
        models <- list(models)
    }

    metrics <- metric_set(precision,accuracy,recall,f_meas,roc_auc,pr_auc)

    model_metrics <- do.call(rbind,lapply(models,FUN = function(model){

        modelClass <- class(extract_spec_parsnip(model))[1]

        metrics <- model %>%
            collect_predictions() %>%
            metrics(truth = use,estimate = .pred_class,.pred_FALSE) %>%
            mutate(model = modelClass)
    }))

    return(model_metrics)
}
