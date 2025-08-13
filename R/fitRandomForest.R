#' fitRandomForest
#'
#' @param data model data
#' @param model model recipe object
#' @param folds crossfold validation object
#' @param metric model performance metric (Default: "accuracy")
#' @param importance Ranger RF engine importance parameter (Default: "impurity")
#' @param trees number of trees (Default: 1000)
#'
#' @return fitted model
#' @export
#'
fitRandomForest <- function(data = NULL,model = NULL,folds = NULL,
                            metric = "accuracy",importance = "impurity",
                            trees = 1000){

    rlang::arg_match(metric,c("precision","accuracy","recall",
                              "f_meas","roc_auc","pr_auc"),multiple = F)

    ## Set randomforest using ranger engine/function
    rf_ranger <- parsnip::rand_forest(mtry = tune::tune(),min_n = tune::tune(),trees = trees) %>%
        parsnip::set_engine("ranger",importance = importance) %>%
        parsnip::set_mode("classification")

    ## Set workflow by combining model and recipe into single model object
    rf_ranger_WF <- workflows::workflow() %>%
        workflows::add_model(rf_ranger) %>%
        workflows::add_recipe(model)



    rf_res <- rf_ranger_WF %>%
        tune::tune_grid(folds,
                  grid = 25,
                  control = tune::control_grid(save_pred = TRUE),
                  metrics = yardstick::metric_set(precision,accuracy,recall,f_meas,roc_auc,pr_auc))

    rf_best <- rf_res %>%
        tune::select_best(metric = metric,n = 10)

    final_rf <- tune::finalize_workflow(rf_ranger_WF,rf_best)
    final_rf_res <- tune::last_fit(final_rf, data)
    return(final_rf_res)
}
