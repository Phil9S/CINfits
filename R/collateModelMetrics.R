collateModelMetrics <- function(models){
    use=`.pred_class`=`.pred_FALSE`=NULL
    if(is.null(models)){
        stop("no models provided")
    }

    if(!inherits(models,what = "list")){
        models <- list(models)
    }

    metrics <- yardstick::metric_set(yardstick::precision,
                                     yardstick::accuracy,
                                     yardstick::recall,
                                     yardstick::f_meas,
                                     yardstick::roc_auc,
                                     yardstick::pr_auc)

    model_metrics <- do.call(rbind,lapply(models,FUN = function(model){

        modelClass <- class(workflowsets::extract_spec_parsnip(model))[1]

        metrics <- model %>%
            tune::collect_predictions() %>%
            metrics(truth = use,estimate = .pred_class,.pred_FALSE) %>%
            dplyr::mutate(model = modelClass)
    }))

    return(model_metrics)
}
