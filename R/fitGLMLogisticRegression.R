#' fitGLMLogisticRegression
#'
#' @param data model data
#' @param model model recipe object
#' @param folds crossfold validation object
#' @param mixture GLM lasso-ridge regression mixture (Default: 1).
#'  Value of 0 is ridge, 1 is lasso, float values between are elastic net.
#' @param metric model performance metric (Default: "accuracy")
#'
#' @return fitted model
#' @export
fitGLMLogisticRegression <- function(data = NULL,model = NULL,folds = NULL,mixture = 1,metric = "accuracy"){

    rlang::arg_match(metric,c("precision","accuracy","recall",
                              "f_meas","roc_auc","pr_auc"),multiple = F)
    ## Set logistic regression with glmnet engine/function
    lgm_glmnet <- parsnip::logistic_reg(penalty = tune::tune(),mixture = mixture) %>%
        parsnip::set_engine("glmnet")

    ## Set workflow by combining model and recipe into single model object
    lgm_glmnet_WF <- workflows::workflow() %>%
        workflows::add_model(lgm_glmnet) %>%
        workflows::add_recipe(model)

    lr_reg_grid <- tidyr::tibble(penalty = 10^seq(-4, -1, length.out = 30))

    lr_res <- lgm_glmnet_WF %>%
        tune::tune_grid(resamples = folds,
                  grid = lr_reg_grid,
                  control = tune::control_grid(save_pred = TRUE),
                  metrics = yardstick::metric_set(yardstick::precision,
                                                  yardstick::accuracy,
                                                  yardstick::recall,
                                                  yardstick::f_meas,
                                                  yardstick::roc_auc,
                                                  yardstick::pr_auc))

    lr_best <- lr_res %>%
        tune::select_best(metric = metric,n = 10)

    lgm_glmnet_WF <- tune::finalize_workflow(lgm_glmnet_WF,lr_best)

    lr_res_final <- tune::last_fit(lgm_glmnet_WF,data)
    return(lr_res_final)
}
