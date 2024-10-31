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
    modelRecipe <- recipes::recipe(use ~ .,data = data) %>%
        recipes::update_role(sample,new_role = "ID") %>%
        recipes::step_zv(recipes::all_predictors()) %>%
        recipes::step_corr(recipes::all_predictors())

    folds <- rsample::vfold_cv(data,
                      v = folds,
                      strata = strata)

    return(list(recipe = modelRecipe,cv = folds))
}
