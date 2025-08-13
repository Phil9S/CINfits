#' fitInteractive
#'
#' Launch an interactive shiny application  to perform copy number fitting and
#' quality control
#'
#' @param data data.frame or list of segmented copy number profiles
#' @param metadata data.frame containing meta data for the samples contained in
#'   data. This should include at least 3 columns; 'sample', 'ploidy', and
#'   'purity'.
#' @param autoSave boolean value to set if QC table should be automatically saved
#'   to a temporary location during interactive fitting.
#' @param autoSaveInt Time increment (minutes) at which to perform an auto save
#'
#' @return interactive shiny application
#' @export
#'
fitInteractive <- function(data=NULL,metadata=NULL,autoSave=FALSE,autoSaveInt=15){
    ## R CMD check notes fix
    predict=specificity=sensitivity=`.estimator`=accuracy=notes=NULL
    notes=model=metric_set=precision=use=savedModelMetrics=NULL

    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop(
            "Package \"shiny\" must be installed to use interactive fitting",
            call. = FALSE
        )
    }
    if (!requireNamespace("DT", quietly = TRUE)) {
        stop(
            "Package \"DT\" must be installed to use interactive fitting",
            call. = FALSE
        )
    }

    if (!requireNamespace("shinyjs", quietly = TRUE)) {
        stop(
            "Package \"shinyjs\" must be installed to use interactive fitting",
            call. = FALSE
        )
    }

    if (!requireNamespace("shinycssloaders", quietly = TRUE)) {
        stop(
            "Package \"shinycssloaders\" must be installed to use interactive fitting",
            call. = FALSE
        )
    }

    if (!requireNamespace("colourpicker", quietly = TRUE)) {
        stop(
            "Package \"colourpicker\" must be installed to use interactive fitting",
            call. = FALSE
        )
    }
    if(is.null(data)){
        stop("no data")
    }
    if(is.null(metadata)){
        stop("no metadata")
    }

    model_choices <- c("x-gradient-boost-trees",
                       "random-forest",
                       "logistic-regression-ridge",
                       "logistic-regression-lasso",
                       "logistic-regression-elasticNet",
                       "support-vector-machine-RBF")

    metric_choices <- c("precision","accuracy","recall",
                        "f_meas","roc_auc","pr_auc")

    ## set seed
    set.seed("0990")

    ## Initial qc table
    qctable <- generateQCTable(data,metadata)

    if(all(c("nAraw","nBraw") %in% colnames(data))){
        AS <- TRUE
    } else {
        AS <- FALSE
    }

    # drop 1bp segments
    data <- data[which(data$end - data$start > 1),]

    data <- droplevels(data)

    data.list <- split(data,f = data$sample)

    ui <- shiny::fluidPage(
        shinyjs::useShinyjs(),
        shiny::tags$head(shiny::tags$style(".btnD { vertical-align: middle; width: 100%; margin: 10px;}")),
        shiny::navbarPage(title ="CNfits",collapsible = TRUE,position = "fixed-top",
            shiny::tabPanel(title = "Fittng",value = "fitTab",
                shiny::fluidRow(
                    shiny::h3(shiny::tags$b("Copy number fitting"))
                ),
                shiny::fluidRow(
                    shiny::column(5,
                        shiny::fluidRow(
                            shiny::textOutput(outputId = "autosave"),
                            shiny::textOutput(outputId = "message"),
                        ),
                        shiny::fluidRow(
                            shiny::column(5,
                                shiny::h4("Sample"),
                                shiny::selectInput("var",label = NULL,choices = NULL),
                                shiny::textOutput("sampleN"),
                                shiny::fluidRow(
                                    shiny::column(4,shiny::checkboxInput("as",label = "allele-specific",value = FALSE)),
                                    shiny::column(4,shiny::checkboxInput("round_values",label = "round segments",value = FALSE)),
                                    shiny::column(4,shiny::checkboxInput("smoothProfile",label = "smooth segments",value = FALSE))
                                ),
                                shiny::conditionalPanel(condition = "input.smoothProfile == true",
                                    shiny::h4("Smoothing parameters"),
                                    shiny::sliderInput("smoothFactor",label = "smoothing factor",min = 0,max = 0.5,step = 0.01,value = 0.12)),
                                shiny::h4("Upload QC"),
                                shiny::fileInput("uploadQC",label = NULL,multiple = F,
                                                 accept = c(".tsv",".csv"),
                                                 buttonLabel = "upload",
                                                 placeholder = "upload qc file...")
                            ),
                            shiny::column(5,offset = 1,
                                shiny::fluidRow(
                                    shiny::h4("Fit parameters"),
                                    shiny::column(6,shiny::sliderInput("pl_new","ploidy",min = 1,max = 8,step = 0.01,value = 2)),
                                    shiny::column(6,shiny::sliderInput("pu_new","purity",min = 0.2,max = 1,step = 0.01,value = 0.7))
                                ),
                                shiny::fluidRow(
                                    shiny::actionButton(inputId = "resetPlPu",label = "reset fit parameters",width = '100%',height = '100px')
                                ),
                                shiny::fluidRow(
                                shiny::h4("Fit selection"),
                                    shiny::fluidRow(
                                        shiny::column(12,
                                                      shiny::column(4,shiny::actionButton("accept_fit",label = "accept",width = '100%')),
                                                      shiny::column(4,shiny::actionButton("reject_fit",label = "reject",width = '100%')),
                                                      shiny::column(4,shiny::actionButton("refit_fit",label = "refit",width = '100%'))
                                        )
                                    ),
                                shiny::br(),
                                    shiny::column(12,
                                                  shiny::textAreaInput("notes_fit",label = NULL,placeholder = "Add fit notes here...",width = '100%')
                                    )
                                )
                            )
                        ),
                        shiny::fluidRow(

                        ),
                        shiny::fluidRow(
                            shiny::column(12,
                                shiny::h4("QC table"),
                                shiny::fluidRow(
                                    DT::dataTableOutput(outputId = "qc_table")
                                ),
                                shiny::fluidRow(
                                    shiny::column(2,offset = 8,shiny::downloadButton(outputId = "saveQC",class = "btnD",label = "save QC")),
                                    shiny::column(2,shiny::downloadButton(outputId = "saveFits",class = "btnD",label = "save fits"))
                                )
                            )
                        )
                    ),
                    shiny::column(6,offset = 1,
                        shiny::fluidRow(
                            shiny::h3("Original profile"),
                                shiny::column(9,shiny::plotOutput("original_fit")),
                                shiny::column(2,shiny::tableOutput("fit_stat"))
                        ),
                        shiny::fluidRow(
                            shiny::h3("refit profile"),
                                shiny::column(9,shiny::plotOutput("new_fit")),
                                shiny::column(2,shiny::tableOutput("new_stat"))
                        )
                    )
                )
            ),
            shiny::tabPanel(title = "ML",value = "mlTab",
                shiny::fluidRow(
                    shiny::h3(shiny::tags$b("Automate fitting (ML)"))
                ),
                shiny::fluidRow(
                    shiny::fluidRow(
                        shiny::column(7,
                            shiny::h4("Model selection"),
                            shiny::column(11,offset = 1,
                                shiny::fluidRow(
                                    shiny::fluidRow(
                                        shiny::fluidRow(
                                            shiny::column(4,
                                                shiny::selectInput("mlmodel",multiple = F,selectize = TRUE,
                                                   label = "model",choices = c("Select model"="",model_choices))),
                                            shiny::column(2,
                                                shiny::selectInput("mlmetrics",multiple = F,selected = "accuracy",selectize = TRUE,
                                                           label = "metric",choices = c("Select metric"="",metric_choices)),
                                                shiny::actionButton(inputId = "runModel",class = "btnD",label = "Run", style = 'margin:0px')
                                            )
                                        ),
                                        shiny::sliderInput(inputId = "proportion",label = "split proportion",
                                                           min = 0.05,max = 0.95,value = 0.7),
                                        shiny::sliderInput(inputId = "folds",label = "CV folds",
                                                           min = 2,max = 10,value = 5)
                                    ),
                                    shiny::fluidRow(
                                        shiny::column(3,shiny::h4("Available data"),shiny::textOutput("sampleFitted")),
                                        shiny::column(4,shiny::h4("Current split"),shiny::textOutput("datasplit"))
                                    )
                                ),
                                shiny::fluidRow(
                                    shiny::fluidRow(shiny::h4("Model statistics"),shiny::textOutput("modelName")),
                                    shiny::column(12,shinycssloaders::withSpinner(DT::DTOutput("metrics"),caption = "running model"))
                                )
                            )
                        ),
                        shiny::column(5,
                            shiny::h4("Model performance"),
                            shiny::fluidRow(
                                shinycssloaders::withSpinner(shiny::plotOutput("roc_plot"),caption = "running model")
                            ),
                            shiny::fluidRow(
                                #shiny::h4("Model prediction"),
                                shiny::column(3,
                                    shiny::selectInput(inputId = "modelSelected",
                                                       label = "model",
                                                       choices = c("Select model"=""),
                                                       multiple = F),
                                    shiny::actionButton(inputId = "applyModel",label = "apply model")#,style = 'margin-top:25px')
                                ),
                                shiny::column(8,
                                    shiny::fluidRow(
                                        DT::DTOutput(outputId = "predictTable")
                                    )
                                )
                            ),
                            shiny::fluidRow(
                                shiny::column(2,offset = 7,
                                    shiny::downloadButton(outputId = "saveModel",label = "save model")
                                ),
                                shiny::column(2,
                                    shiny::downloadButton(outputId = "savePredict",label = "save table")
                                )
                            )
                        )
                    )
                )
            ),
            shiny::tabPanel(title = "options",value = "optsTab",
                shiny::fluidRow(
                    shiny::h3(shiny::tags$b("Options"))
                ),
                shiny::fluidRow(
                    shiny::h4("Segment colours"),
                    shiny::conditionalPanel(condition = "input.as == true",
                        shiny::column(6,colourpicker::colourInput(inputId = "nA",
                                                                  label = "allele A",
                                                                  value = "blue",allowTransparent = F,
                                                                  showColour = "background",
                                                                  closeOnClick = TRUE)),
                        shiny::column(6,colourpicker::colourInput(inputId = "nB",
                                                                  label = "allele B",
                                                                  value = "darkred",
                                                                  allowTransparent = F,
                                                                  showColour = "background",
                                                                  closeOnClick = TRUE))
                    ),
                    shiny::conditionalPanel(condition = "input.as == false",
                        shiny::column(12,colourpicker::colourInput(inputId = "totcol",
                                                                   label = "total copy number",
                                                                   value = "blue",
                                                                   allowTransparent = F,
                                                                   showColour = "background",
                                                                   closeOnClick = TRUE))
                    )
                )
            )
        )
    )

    server <- function(input, output, session) {

        # disable allele-specific where only total available
        shiny::observe(if(AS == FALSE){
            shinyjs::disable(id = "as")
        })


        shiny::observe(if(input$modelSelected == ""){
            shinyjs::disable(id = "saveModel")
            shinyjs::disable(id = "savePredict")
        } else {
            shinyjs::enable(id = "saveModel")
            shinyjs::enable(id = "savePredict")
        })


        # switch between allele-specific and total CN data
        shiny::observe(if(AS == TRUE){
            shiny::updateCheckboxInput(inputId = "as",value = TRUE)
        })

        # observe QC table
        qcData <- shiny::reactiveValues(data=NULL)
        shiny::observe({
            qcData$data <- qctable
        })

        # Initialise sample list
        shiny::observe({
            shiny::updateSelectizeInput(session,inputId = "var",server = TRUE,
                                     choices = names(data.list))
            output$sampleN <- shiny::renderText({
                currn <- which(names(data.list) %in% input$var)
                messn <- paste0(currn," of ",length(names(data.list)))
                return(messn)
            })
        })
        # upload partial or pre-existing qc file
        shiny::observeEvent(input$uploadQC,{
            inFile <- input$uploadQC
            qcUploadTab <- data.table::fread(inFile$datapath)
            if(!all(colnames(qcUploadTab) %in% c("sample","segments","clonality",
                                                "ploidy","purity","homozygousLoss",
                                                "use","notes"))){
                #output$message <- shiny::renderText("incompatible QC file")
                shiny::showNotification("incompatible QC file",
                                        duration = 5,type = "error",closeButton = T)
            } else {
                qcData$data <- qcUploadTab

                if(any(is.na(qcData$data$use))){
                    min_pos <- min(which(is.na(qcData$data$use)))
                } else {
                    min_pos <- 1
                }
                selectedSample <- names(data.list)[min_pos]
                shiny::updateSelectInput(session,inputId = "var",
                                         selected = selectedSample)
            }

        })

        shiny::observe({
            output$sampleFitted <- shiny::renderText({
                currn <- sum(!is.na(qcData$data$use))
                messn <- paste0("Total fitted: ",currn," of ",length(qcData$data$use))
                return(messn)
            })
        })

        shiny::observe({
            output$datasplit <- shiny::renderText({
                currn <- sum(!is.na(qcData$data$use))
                if(currn < 3){
                    messn <- "insufficent samples"
                    return(messn)
                } else {
                    df <- qcData$data
                    modelData <- df[!is.na(df$use),]
                    datasplit <- splitFittingData(data = modelData,prop = input$proportion)
                    tr <- nrow(datasplit[[2]])
                    te <- nrow(datasplit[[3]])
                    messn <- paste0("training: ",tr," | testing: ",te," | total: ",sum(tr,te))
                }
                return(messn)
            })
        })

        modelMetrics <- shiny::eventReactive(input$runModel,{

            df <- qcData$data
            dfFitted <- df %>%
                dplyr::select(-notes) %>%
                tidyr::drop_na() %>%
                dplyr::mutate(use = factor(use))

            shiny::req(input$mlmodel,nrow(dfFitted)*input$proportion > input$folds)

            dfSplit <- splitFittingData(data = dfFitted,
                                        prop = input$proportion,
                                        strata = "use")

            dsplit <- dfSplit$dataSplit
            ## Set training and test data sets
            trainingData <- dfSplit$trainingData
            testingData <- dfSplit$testingDat

            ModelRecipe_CV <- makeModelRecipe(data = trainingData,folds = input$folds)
            modelRecipe <- ModelRecipe_CV$recipe

            folds <- ModelRecipe_CV$cv
            metric <- input$mlmetrics
            model_name <- input$mlmodel

            switch(input$mlmodel,
                   "x-gradient-boost-trees"={
                       fittedModel <- fitXGBtree(data = dsplit,
                                                 model = modelRecipe,
                                                 folds = folds,
                                                 metric = metric,
                                                 trees = 1000)
                   },
                   "random-forest"={
                        fittedModel <- fitRandomForest(data = dsplit,
                                                       model = modelRecipe,
                                                       folds = folds,
                                                       metric = metric,
                                                       importance = "impurity",
                                                       trees = 1000)
                   },
                   "logistic-regression-ridge"={
                       fittedModel <- fitGLMLogisticRegression(data = dsplit,
                                                               model = modelRecipe,
                                                               folds = folds,
                                                               mixture = 0,
                                                               metric = metric)
                   },
                   "logistic-regression-lasso"={
                       fittedModel <- fitGLMLogisticRegression(data = dsplit,
                                                               model = modelRecipe,
                                                               folds = folds,
                                                               mixture = 1,
                                                               metric = metric)
                   },
                   "logistic-regression-elasticNet"={
                        fittedModel <- fitGLMLogisticRegression(data = dsplit,
                                                                model = modelRecipe,
                                                                folds = folds,
                                                                mixture = 0.5,
                                                                metric = metric)
                   },
                   "support-vector-machine-RBF"={
                        fittedModel <-  fitSVM(data = dsplit,
                                               model = modelRecipe,
                                               folds = folds,
                                               metric = metric)
                   }
            )

            metrics <- data.frame(metric=metric,split=input$proportion,folds=input$folds)
            if(is.null(savedModelMetrics)){
                #print("new")
                savedModelMetrics <<- list(list(model=fittedModel,model_name=model_name,metric=metrics))
            } else {
                #print("append")
                newModel <- list(model=fittedModel,model_name=model_name,metric=metrics)
                savedModelMetrics <<- append(savedModelMetrics,list(newModel))
            }

            return(savedModelMetrics)
        })

        shiny::observe({
            shiny::req(input$runModel)
            models <- metricsTable()
            shiny::updateSelectInput(inputId = "modelSelected",
                                     choices = models$model,
                                     selected = NULL)
        })

        metricsTable <- shiny::reactive({
            modelData <- modelMetrics()
            #print(length(modelData))

            fittedModels <- lapply(modelData,FUN = function(x){
                x[["model"]]
            })

            model_name <- do.call(rbind,lapply(1:length(modelData),FUN = function(x){
                y <- modelData[[x]]
                y <- y[["model_name"]]
                #print(y)
                z <- data.frame(model=rep(paste0("model_",x),length(metric_choices)),
                                model_type=rep(y,length(metric_choices)))
                #print(z)
                return(z)
            }))

            # model_name <- modelMetrics()[["model_name"]]
            metric_name <- do.call(rbind,lapply(modelData,FUN = function(x){
                y <- x[["metric"]]
                y <- y[rep(1,length(metric_choices)),]
                return(y)
            }))

            #print(model_name)
            #print(metric_name)
            tb <- collateModelMetrics(fittedModels)
            tb$model <- model_name$model
            tb$model_type <- model_name$model_type
            tb <- cbind(tb,metric_name)
            #print(tb)

            tb <- tb %>%
                dplyr::select(-.estimator) %>%
                tidyr::pivot_wider(names_from = ".metric",
                                   id_cols = c("model","model_type","metric","folds","split"),
                                   values_from = ".estimate") %>%
                dplyr::relocate(model,.before = 1) %>%
                dplyr::mutate(dplyr::across(dplyr::where(is.numeric),~round(.x,digits = 2)))

            return(tb)
        })

        output$metrics <- DT::renderDT({

            DT::datatable(metricsTable(),rownames = FALSE,extensions = 'Scroller',
                          options = list(
                              class = "display",
                              pageLength = 5,
                              ordering=T,
                              searching=F,
                              scrollY = 300,
                              scroller = TRUE,
                              searchHighlight = T))
        })

        rocMetrics <- shiny::reactive({
            modelData <- modelMetrics()

            fittedModels <- lapply(modelData,FUN = function(x){
                x[["model"]]
            })

            model_name <- do.call(rbind,lapply(1:length(modelData),FUN = function(x){
                m <- modelData[[x]]
                y <- m[["model_name"]]
                #print(m)
                lenrep <- length(unique(m$model$.predictions[[1]]$.pred_FALSE)) + 2
                #print(lenrep)
                z <- data.frame(model=rep(paste0("model_",x),times=lenrep),
                                model_type=rep(y,times=lenrep))
                return(z)
            }))
            #print(model_name)
            modelROC <- collateModelRoc(fittedModels,usePR = F)
            #print(modelROC)
            modelROC$model <- model_name$model
            modelROC$model_type <- model_name$model_type

            modelROC$model <- factor(modelROC$model,
                                     levels = stringr::str_sort(unique(modelROC$model),
                                                       numeric = T))
            return(modelROC)
        })

        output$roc_plot <- shiny::renderPlot({

            rocplot <- ggplot2::ggplot(rocMetrics(),
                              ggplot2::aes(x = 1 - specificity,
                                  y = sensitivity,
                                  col = model)) +
                ggplot2::geom_path(lwd = 1.5, alpha = 0.8) +
                ggplot2::geom_abline(lty = 3) +
                ggplot2::coord_equal() +
                ggplot2::theme_bw()

            return(rocplot)

        })

        predictOut <- shiny::eventReactive(input$applyModel,{
            qcTablePred <- qcData$data
            #qcTablePred <- qcTablePred[is.na(qcTablePred$use),]

            modelIndex <- as.numeric(gsub(pattern = "model_",replacement = "",x = input$modelSelected))
            #print(modelIndex)
            selectedModelList <- modelMetrics()[[modelIndex]]
            selectedModel <- selectedModelList[["model"]]
            #print(selectedModel)

            predModel <- workflowsets::extract_workflow(selectedModel)
            #print(predModel)
            #print(qcTablePred)
            newClass <- stats::predict(predModel,qcTablePred,type = "class")
            newProb <- round(stats::predict(predModel,qcTablePred,type = "prob"),digits = 2)

            qcTablePred <- qcTablePred[,c("sample","use")]
            qcTablePred <- cbind(qcTablePred,newClass,newProb)

            return(qcTablePred)
        })

        output$predictTable <- DT::renderDT({
            DT::datatable(predictOut(),rownames = FALSE,extensions = 'Scroller',
                          options = list(
                              class = "display",
                              pageLength = 5,
                              ordering=F,
                              searching=F,
                              scrollY = 300,
                              scroller = TRUE,
                              searchHighlight = T))
        })


        # action on accepting fit
        shiny::observeEvent(input$accept_fit,{

            df <- qcData$data
            df$use[df$sample == input$var] <- TRUE

            df$notes[df$sample == input$var] <- input$notes_fit
            shiny::updateTextAreaInput(inputId = "notes_fit",value = "")

            qcData$data <- df

            if(any(is.na(qcData$data$use))){
                min_pos <- min(which(is.na(qcData$data$use)))
            } else {
                min_pos <- 1
            }
            selectedSample <- names(data.list)[min_pos]

            shiny::updateSelectInput(session,inputId = "var",
                                     selected = selectedSample)
        })

        # action on accepting refit
        shiny::observeEvent(input$refit_fit,{
            df <- qcData$data
            df$use[df$sample == input$var] <- TRUE

            orig_ploidy <- qcData$data$ploidy[qcData$data$sample == input$var]
            orig_purity <- qcData$data$purity[qcData$data$sample == input$var]

            df$notes[df$sample == input$var] <- paste0("REFIT: original ploidy=",
                                                       orig_ploidy,
                                                       ", original purity=",
                                                       orig_purity,". ",input$notes_fit)

            df$ploidy[df$sample == input$var] <- input$pl_new
            df$purity[df$sample == input$var] <- input$pu_new
            shiny::updateTextAreaInput(inputId = "notes_fit",value = "")

            ## update profile in data.list
            new_fit <- rescaleFit(data = data.list[[input$var]],
                                  old_ploidy = orig_ploidy,
                                  old_purity = orig_purity,
                                  new_ploidy = input$pl_new,
                                  new_purity = input$pu_new,
                                  alleleSpecific=input$as)

            data.list[[input$var]] <<- new_fit

            qcData$data <- df

            if(any(is.na(qcData$data$use))){
                min_pos <- min(which(is.na(qcData$data$use)))
            } else {
                min_pos <- 1
            }
            selectedSample <- names(data.list)[min_pos]
            shiny::updateSelectInput(session,inputId = "var",
                                     selected = selectedSample)

        })

        # action on reject fit
        shiny::observeEvent(input$reject_fit,{

            df <- qcData$data
            df$use[df$sample == input$var] <- FALSE

            df$notes[df$sample == input$var] <- input$notes_fit
            shiny::updateTextAreaInput(session,inputId = "notes_fit",value = "")

            qcData$data <- df

            if(any(is.na(qcData$data$use))){
                min_pos <- min(which(is.na(qcData$data$use)))
            } else {
                min_pos <- 1
            }
            selectedSample <- names(data.list)[min_pos]

            shiny::updateSelectInput(session,inputId = "var",
                                     selected = selectedSample)
        })

        # Update sliders to provided fit
        shiny::observe({
            shiny::updateSliderInput(session,inputId = "pl_new",
                              value = qcData$data$ploidy[qcData$data$sample == input$var])
            shiny::updateSliderInput(session,inputId = "pu_new",
                              value = qcData$data$purity[qcData$data$sample == input$var])
        })

        ## Allow for resetting fitting values back to original values
        shiny::observeEvent(input$resetPlPu,{
            shiny::updateSliderInput(session,inputId = "pl_new",
                                     value = qcData$data$ploidy[qcData$data$sample == input$var])
            shiny::updateSliderInput(session,inputId = "pu_new",
                                     value = qcData$data$purity[qcData$data$sample == input$var])
        })

        output$original_fit <- shiny::renderPlot({
            shiny::req(input$var)
            orig_fit <- data.list[[input$var]]
            orig_purity <- qcData$data$purity[qcData$data$sample == input$var]

            if(input$round_values){
                orig_fit$segVal <- round(orig_fit$segVal)
                if(input$as){
                    orig_fit$nAraw <- round(orig_fit$nAraw)
                    orig_fit$nBraw <- round(orig_fit$nBraw)
                }
            }

            if(input$smoothProfile){
                orig_fit <- smoothProfile(orig_fit,
                                        smoothingFactor = input$smoothFactor,
                                        alleleSpecific = input$as)
            }

            if(input$as){
                cols <- c(input$nA,input$nB)
            } else {
                cols <- input$totcol
            }

            plotProfile(orig_fit,sample = input$var,
                        cn.max = 15,purity = orig_purity,
                        alleleSpecific = input$as,cols = cols)
        })

        output$new_fit <- shiny::renderPlot({
            shiny::req(input$var)
            orig_ploidy <- qcData$data$ploidy[qcData$data$sample == input$var]
            orig_purity <- qcData$data$purity[qcData$data$sample == input$var]

            new_fit <- rescaleFit(data = data.list[[input$var]],
                                  old_ploidy = orig_ploidy,
                                  old_purity = orig_purity,
                                  new_ploidy = input$pl_new,
                                  new_purity = input$pu_new,
                                  alleleSpecific=input$as)

            if(input$round_values){
                new_fit$segVal <- round(new_fit$segVal)
                if(input$as){
                    new_fit$nAraw <- round(new_fit$nAraw)
                    new_fit$nBraw <- round(new_fit$nBraw)
                }
            }

            if(input$smoothProfile){
                new_fit <- smoothProfile(new_fit,
                                         smoothingFactor = input$smoothFactor,
                                         alleleSpecific = input$as)
            }

            if(input$as){
                cols <- c(input$nA,input$nB)
            } else {
                cols <- input$totcol
            }

            plotProfile(new_fit,sample = input$var,
                        cn.max = 15,purity = input$pu_new,
                        alleleSpecific = input$as,cols = cols)
        })

        output$fit_stat <- shiny::renderTable({
            fitTab <- calculateCINStats(data.list)
            fitTab <- as.data.frame(fitTab[rownames(fitTab) == input$var,])
            fitTab <- data.frame(statistic=rownames(fitTab),value=fitTab[,1])
            return(fitTab)
        })

        output$new_stat <- shiny::renderTable({
            shiny::req(input$var)
            orig_ploidy <- qcData$data$ploidy[qcData$data$sample == input$var]
            orig_purity <- qcData$data$purity[qcData$data$sample == input$var]

            newTab <- rescaleFit(data = data.list[[input$var]],
                                 old_ploidy = orig_ploidy,
                                 old_purity = orig_purity,
                                 new_ploidy = input$pl_new,
                                 new_purity = input$pu_new,
                                 alleleSpecific=input$as)

            if(input$smoothProfile){
                newTab <- smoothProfile(newTab,
                                        smoothingFactor = input$smoothFactor,
                                        alleleSpecific = input$as)
            }

            newTab <- calculateCINStats(newTab)
            newTab <- as.data.frame(newTab[rownames(newTab) == input$var,])
            newTab <- data.frame(statistic=rownames(newTab),value=newTab[,1])
            return(newTab)
        })

        output$qc_table <- DT::renderDataTable({
            DT::datatable(qcData$data,rownames = FALSE,extensions = 'Scroller',
                      options = list(
                          columnDefs = list(list(targets = "_all",render = DT::JS(
                              "function(data, type, row, meta) {",
                              "return type === 'display' && data != null && data.length > 11 ?",
                              "'<span title=\"' + data + '\">' + data.substr(0, 7) + '...</span>' : data;",
                              "}")
                                    )),
                          class = "display",
                          pageLength = 5,
                          ordering=F,
                          scrollY = 200,
                          scroller = TRUE,
                          searchHighlight = T
                                     ))
        })

        if(autoSave){
            autoCount <- 0
            autoSaveTrigger <- shiny::reactiveTimer(intervalMs = autoSaveInt*60000,session = session)
            shiny::observe({
                autoSaveTrigger()
                if(autoCount == 0){
                    output$autosave <- shiny::renderText({paste0("Last saved : never")})
                    autoCount <- 1
                }
                wd <- getwd()
                time <- gsub(":","-",gsub(" ","_",as.character(Sys.time())))

                # save QC
                file <- paste0(wd,"/autosave_",time,"_QC_table",".tsv")
                utils::write.table(shiny::isolate(qcData$data),file,quote = F,sep = "\t",append = F,
                            row.names = F,col.names = T)

                # Save fits
                dataF <- do.call(rbind,lapply(data.list,FUN = function(x){
                    return(x)}
                ))
                fileF <- paste0(wd,"/autosave_",time,"_seg_table",".tsv")
                utils::write.table(dataF,fileF,quote = F,sep = "\t",append = F,
                            row.names = F,col.names = T)

                # report last save
                lastTime <- round(Sys.time())
                output$autosave <- shiny::renderText({paste0("Last saved : ",as.character(lastTime))})
            })
        }

        output$saveQC <- shiny::downloadHandler(
            filename = function() {
                paste0(Sys.Date(),"_QC_table",".tsv")
            },
            content = function(file) {
                utils::write.table(qcData$data,file,quote = F,sep = "\t",append = F,
                            row.names = F,col.names = T)
            }
        )

        output$saveModel <- shiny::downloadHandler(
            filename = function() {
                modelIndex <- as.numeric(gsub(pattern = "model_",
                                              replacement = "",
                                              x = input$modelSelected))
                selectedModelList <- modelMetrics()[[modelIndex]]
                modelName <- selectedModelList[["model_name"]]
                paste0(Sys.Date(),"_CN_fitting_model_",modelName,".rds")
            },
            content = function(file) {
                modelIndex <- as.numeric(gsub(pattern = "model_",
                                              replacement = "",
                                              x = input$modelSelected))
                #print(modelIndex)
                selectedModelList <- modelMetrics()[[modelIndex]]
                selectedModel <- selectedModelList[["model"]]
                #print(selectedModel)
                predModel <- workflowsets::extract_workflow(selectedModel)
                selectedModelList$extracted_model <- predModel

                saveRDS(selectedModelList,file)
            }

        )

        output$savePredict <- shiny::downloadHandler(
            filename = function() {
                modelIndex <- as.numeric(gsub(pattern = "model_",
                                              replacement = "",
                                              x = input$modelSelected))
                selectedModelList <- modelMetrics()[[modelIndex]]
                modelName <- selectedModelList[["model_name"]]
                paste0(Sys.Date(),"_QC_table_fitted_",modelName,".tsv")
            },
            content = function(file) {
                predTable <- predictOut()

                predTable <- dplyr::left_join(qcData$data,predTable,by = c("sample","use"))

                utils::write.table(predTable,file,quote = F,sep = "\t",append = F,
                                   row.names = F,col.names = T)
            }

        )

        output$saveFits <- shiny::downloadHandler(
            filename = function() {
                paste0(Sys.Date(),"_seg_table",".tsv")
            },
            content = function(file) {
                dataF <- do.call(rbind,lapply(data.list,FUN = function(x){
                    return(x)}
                    ))
                utils::write.table(dataF,file,quote = F,sep = "\t",append = F,
                            row.names = F,col.names = T)
            }
        )
    }
    shiny::shinyApp(ui, server,options = list(launch.browser = TRUE))
}

# END
