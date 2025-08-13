#' fitSVM
#'
#' @param data model data
#' @param model model recipe object
#' @param folds crossfold validation object
#' @param metric model performance metric (Default: "accuracy")
#'
#' @return fitted model
#' @export
fitSVM <- function(data = NULL,model = NULL,folds = NULL,metric = "accuracy"){

    rlang::arg_match(metric,c("precision","accuracy","recall",
                              "f_meas","roc_auc","pr_auc"),multiple = F)

    #set RBF svm
    linsvm_mod <- parsnip::svm_rbf(cost = tune::tune(), rbf_sigma = tune::tune()) %>%
        parsnip::set_mode("classification") %>%
        parsnip::set_engine("kernlab")

    ## Set workflow by combining model and recipe into single model object
    linsvm_WF <- workflows::workflow() %>%
        workflows::add_model(linsvm_mod) %>%
        workflows::add_recipe(model)

    #lr_reg_grid <- tidyr::tibble(penalty = 10^seq(-4, -1, length.out = 30))
    linsvm_res <- linsvm_WF %>%
        tune::tune_grid(resamples = folds,
                        #grid = lr_reg_grid,
                        control = tune::control_grid(save_pred = TRUE),
                        metrics = yardstick::metric_set(yardstick::precision,
                                                        yardstick::accuracy,
                                                        yardstick::recall,
                                                        yardstick::f_meas,
                                                        yardstick::roc_auc,
                                                        yardstick::pr_auc))

    linsvm_best <- linsvm_res %>%
        tune::select_best(metric,n = 10)

    linsvm_WF <- tune::finalize_workflow(linsvm_WF,linsvm_best)

    linsvm_res_final <- tune::last_fit(linsvm_WF,data)
    return(linsvm_res_final)
}
