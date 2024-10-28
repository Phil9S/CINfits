#' fitSVM
#'
#' @param model model recipe object
#' @param folds crossfold validation object
#' @param mixture GLM lasso-ridge regression mixture (Default: 1).
#'  Value of 0 is ridge, 1 is lasso, float values between are elastic net.
#' @param metric model performance metric
#'
#' @return fitted model
#' @export
fitSVM <- function(data = NULL,model = NULL,folds = NULL,mixture = 1,metric = "accuracy"){
    #set linear svm
    linsvm_mod <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
        set_mode("classification") %>%
        set_engine("kernlab")

    ## Set workflow by combining model and recipe into single model object
    linsvm_WF <- workflow() %>%
        add_model(linsvm_mod) %>%
        add_recipe(model)

    #lr_reg_grid <- tidyr::tibble(penalty = 10^seq(-4, -1, length.out = 30))

    lr_res <- lgm_glmnet_WF %>%
        tune::tune_grid(resamples = folds,
                        #grid = lr_reg_grid,
                        control = tune::control_grid(save_pred = TRUE),
                        metrics = yardstick::metric_set(metric))

    linsvm_best <- linsvm_res %>%
        select_best(metric,n = 10)

    linsvm_WF <- finalize_workflow(linsvm_WF,linsvm_best)

    linsvm_res_final <- last_fit(linsvm_WF,data)
    return(linsvm_res_final)
}
