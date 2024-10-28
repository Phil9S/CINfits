#' makeModelRecipe
#'
#' @param data training data
#'
#' @return recipe object for tidymodels
#' @export
#'
makeModelRecipe <- function(data = NULL,folds = 10,strata = "use"){
    ## Set recipe for model
    # update sample column to be ID rather than predictor
    modelRecipe <- recipe(use ~ .,data = data) %>%
        update_role(sample,new_role = "ID") %>%
        step_zv(all_predictors()) %>%
        step_corr(all_predictors())

    folds <- vfold_cv(data,
                      v = folds,
                      strata = strata)

    return(list(recipe = modelRecipe,cv = folds))
}
