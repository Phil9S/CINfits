#' splitFittingData
#' Split manually fitted QC data into training and test data
#' @param data QC data table
#' @param prop proportion to retain for training
#' @param strata variable to stratify in order to maintain class balance
#'
#' @return list of sample split, training split, and testing split samples
#' @export
#'
splitFittingData <- function(data = NULL,prop = 0.7,strata = "use"){

    dataSplit <- rsample::initial_split(data,
                                        prop = prop,
                                        strata = strata)
    ## Set training and test data sets
    trainingData <- rsample::training(dataSplit)
    #validationData <- rsample::validation(dataSplit)
    testingData <- rsample::testing(dataSplit)

    return(list(dataSplit = dataSplit,
                trainingData = trainingData,
                testingData = testingData))
}
