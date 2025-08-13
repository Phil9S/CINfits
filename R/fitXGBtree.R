#' fitXGBtree
#'
#' @param data model data
#' @param model model recipe object
#' @param folds crossfold validation object
#' @param metric model performance metric (Default: "accuracy")
#' @param trees number of trees (Default: 1000)
#'
#' @return fitted model
#' @export
#'
fitXGBtree <- function(data = NULL,model = NULL,folds = NULL,
                            metric = "accuracy",trees = 1000){

    rlang::arg_match(metric,c("precision","accuracy","recall",
                              "f_meas","roc_auc","pr_auc"),multiple = F)

    ## Set xbg trees using xgboost engine/function
    bt_xgb <- parsnip::boost_tree(mtry = tune::tune(),
                         tree_depth = tune::tune(),
                         learn_rate = tune::tune(),
                         sample_size = tune::tune(),
                         loss_reduction = tune::tune(),
                         min_n = tune::tune(),
                         trees = trees) %>%
        parsnip::set_mode("classification") %>%
        parsnip::set_engine("xgboost")

    ## Set workflow by combining model and recipe into single model object
    bt_WF <- workflows::workflow() %>%
        workflows::add_model(bt_xgb) %>%
        workflows::add_recipe(model)

    xgb_grid <- dials::grid_latin_hypercube(
        dials::tree_depth(),
        dials::min_n(),
        dials::loss_reduction(),
        sample_size = dials::sample_prop(),
        dials::finalize(dials::mtry(), rsample::training(data)),
        dials::learn_rate(),
        size = 30
    )

    bt_res <- bt_WF %>%
        tune::tune_grid(resamples = folds,
                  grid = xgb_grid,
                  control = tune::control_grid(save_pred = TRUE),
                  metrics = yardstick::metric_set(yardstick::precision,
                                                  yardstick::accuracy,
                                                  yardstick::recall,
                                                  yardstick::f_meas,
                                                  yardstick::roc_auc,
                                                  yardstick::pr_auc))

    bt_best <- bt_res %>%
        tune::select_best(metric = metric,n = 10)

    final_xgb <- tune::finalize_workflow(bt_WF,bt_best)
    final_xgb_res <- tune::last_fit(final_xgb, data)
    return(final_xgb_res)
}
