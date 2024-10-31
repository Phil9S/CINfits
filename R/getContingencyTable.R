#' getContingencyTable
#'
#' @param model
#' @param var
#'
#' @return contingency table
#' @export
#'
getContingencyTable <- function(model = NULL,var = "use"){
    if(!inherits(lr_res_final,c("last_fit"))){
        stop("model required to be last fit")
    }
    variable <- var

    model %>%
        workflowsets::collect_predictions() %>%
        dplyr::count(.pred_class,(!!sym(variable))) %>%
        tidyr::pivot_wider(id_cols = 1,names_from = 2,values_from = n) %>%
        tibble::column_to_rownames(var = ".pred_class")
}
