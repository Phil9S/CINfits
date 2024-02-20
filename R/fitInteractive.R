#' fitInteractive
#'
#' Launch an interactive shiny application  to perform copy number fitting and
#' quality control
#'
#' @param data data.frame or list of segmented copy number profiles
#' @param metadata data.frame containing meta data for the samples contained in
#'   data. This should include at least 3 columns; 'sample', 'ploidy', and
#'   'purity'.
#'
#' @return interactive shiny application
#' @export
#'
fitInteractive <- function(data=NULL,metadata=NULL){
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

    if (!requireNamespace("colourpicker", quietly = TRUE)) {
        stop(
            "Package \"shinyjs\" must be installed to use interactive fitting",
            call. = FALSE
        )
    }
    if(is.null(data)){
        stop("no data")
    }
    if(is.null(metadata)){
        stop("no metadata")
    }

    ## Initial qc table
    qctable <- generateQCTable(data,metadata)

    if(all(c("nAraw","nBraw") %in% colnames(data))){
        AS <- TRUE
    } else {
        AS <- FALSE
    }

    data.list <- split(data,f = data$sample)
    ui <- shiny::fluidPage(
        shinyjs::useShinyjs(),
        title = "CNfits",
        shiny::titlePanel(title = "CNfits"),
        shiny::column(width = 5,
            shiny::fluidRow(
                shiny::column(width = 4,
                    shiny::fluidRow(
                        shiny::h3("Sample"),
                        shiny::selectInput("var",label = NULL,choices = names(data.list)),
                        shiny::fluidRow(
                            shiny::column(4,shiny::checkboxInput("as",label = "allele-specific",value = FALSE)),
                            shiny::column(4,shiny::checkboxInput("round_values",
                                                                 label = "round segments",value = FALSE)),
                            shiny::column(4,shiny::checkboxInput("smoothProfile",
                                                                 label = "smooth segments",value = FALSE))
                            ),
                        shiny::hr()
                    ),
                    shiny::fluidRow(
                        shiny::conditionalPanel(condition = "input.smoothProfile == true",
                                                shiny::h4("smoothing parameters"),
                                                shiny::sliderInput("smoothFactor",label = "smoothing factor",
                                                                   min = 0,max = 0.5,step = 0.01,value = 0.12),
                                                shiny::hr()
                        )
                    ),
                    shiny::fluidRow(
                            shiny::h4("segment colours"),
                            shiny::conditionalPanel(condition = "input.as == true",
                                            shiny::column(6,
                                            colourpicker::colourInput(inputId = "nA",
                                                                      label = "allele A",
                                                                      value = "red",
                                                                      allowTransparent = F)),
                                            shiny::column(6,
                                            colourpicker::colourInput(inputId = "nB",
                                                                      label = "allele B",
                                                                      value = "blue",
                                                                      allowTransparent = F))
                                            ),
                            shiny::conditionalPanel(condition = "input.as == false",
                                                shiny::column(12,
                                                              colourpicker::colourInput(inputId = "totcol",
                                                                                        label = "total CN",
                                                                                        value = "red",
                                                                                        allowTransparent = F))
                                                )
                            ),
                    shiny::hr(),
                    shiny::fluidRow(
                        shiny::fileInput("uploadQC",
                                         label = "Upload QC file",
                                         multiple = F,
                                         accept = c(".tsv",".csv"),
                                         buttonLabel = "upload",
                                         placeholder = "upload qc file..."),
                        shiny::hr(),
                        shiny::h4("fit selection"),
                        shiny::fluidRow(
                            shiny::column(4,shiny::actionButton("accept_fit",label = "accept",width = '100%')),
                            shiny::column(4,shiny::actionButton("refit_fit",label = "refit",width = '100%')),
                            shiny::column(4,shiny::actionButton("reject_fit",label = "reject",width = '100%'))
                        ),
                        shiny::textAreaInput("notes_fit",label = "notes",
                                             placeholder = "Add fit notes here..."),
                        shiny::hr()
                    )
                ),
                shiny::column(width = 7,offset = 1,
                    shiny::sliderInput("pl_new","ploidy",min = 1,max = 8,
                                       step = 0.01,value = 2),
                    shiny::sliderInput("pu_new","purity",min = 0.2,max = 1,
                                       step = 0.01,value = 0.7),
                    shiny::plotOutput("fit_sunrise",click = "sunrise_click"),
                    shiny::tableOutput(outputId = "info")
                ),
            ),
            shiny::fluidRow(
                shiny::column(width = 12,
                    DT::dataTableOutput(outputId = "qc_table"),
                    shiny::downloadButton(outputId = "saveQC",label = "save")
                )
            ),
        ),
        shiny::column(width = 7,
            shiny::fluidRow(
                    shiny::fluidRow(
                        shiny::h3("Original profile"),
                        shiny::column(width = 8,
                            shiny::plotOutput("original_fit")),
                        shiny::column(width = 1,
                            shiny::tableOutput("fit_stat"))
                    ),
                    shiny::fluidRow(
                        shiny::h3("refit profile"),
                        shiny::column(width = 8,
                            shiny::plotOutput("new_fit")),
                        shiny::column(width = 1,
                            shiny::tableOutput("new_stat"))
                    )
                )
        )
    )

    server <- function(input, output, session) {

        # disable allele-specific where only total available
        shiny::observe(if(AS == FALSE){
            shinyjs::disable(id = "as")
        })

        shiny::observe(if(AS == TRUE){
            shiny::updateCheckboxInput(inputId = "as",value = TRUE)
        })

        # observe QC table
        qcData <- shiny::reactiveValues(data=NULL)
        shiny::observe({
            qcData$data <- qctable
        })

        # upload partial or pre-existing qc file
        observeEvent(input$uploadQC,{
            inFile <- input$uploadQC
            qcUploadTab <- data.table::fread(inFile$datapath)
            qcData$data <- qcUploadTab

            min_pos <- min(which(is.na(qcData$data$use)))
            selectedSample <- names(data.list)[min_pos]
            shiny::updateSelectInput(session,inputId = "var",
                                     selected = selectedSample)
        })

        # action on accepting fit
        shiny::observeEvent(input$accept_fit,{
            #samplePos(samplePos()+1)
            #selectedSample <- names(data.list)[samplePos()]

            df <- qcData$data
            df$use[df$sample == input$var] <- TRUE

            df$notes[df$sample == input$var] <- input$notes_fit
            shiny::updateTextAreaInput(inputId = "notes_fit",value = "")

            qcData$data <- df

            min_pos <- min(which(is.na(qcData$data$use)))
            #if(min_pos < samplePos()){
            #samplePos(min_pos)
            #}
            selectedSample <- names(data.list)[min_pos]

            shiny::updateSelectInput(session,inputId = "var",
                                     selected = selectedSample)
        })

        # action on accepting refit
        shiny::observeEvent(input$refit_fit,{
            # samplePos(samplePos()+1)
            # selectedSample <- names(data.list)[samplePos()]
            #
            # df <- qcData$data
            # df$use[df$sample == input$var] <- TRUE
            #
            # df$notes[df$sample == input$var] <- input$notes_fit
            # updateTextAreaInput(inputId = "notes_fit",value = "")
            #
            # qcData$data <- df
            #
            # shiny::updateSelectInput(inputId = "var",selected = selectedSample)

        })

        # action on reject fit
        shiny::observeEvent(input$reject_fit,{

            #samplePos(samplePos()+1)

            df <- qcData$data
            df$use[df$sample == input$var] <- FALSE

            df$notes[df$sample == input$var] <- input$notes_fit
            shiny::updateTextAreaInput(session,inputId = "notes_fit",value = "")

            qcData$data <- df

            min_pos <- min(which(is.na(qcData$data$use)))
            #if(min_pos < samplePos()){
            #samplePos(min_pos)
            #}
            selectedSample <- names(data.list)[min_pos]

            shiny::updateSelectInput(session,inputId = "var",
                                     selected = selectedSample)
        })

        # Update sliders to provided fit
        shiny::observe({
            shiny::updateSliderInput(session,inputId = "pl_new",
                              value = metadata$ploidy[metadata$sample == input$var])
            shiny::updateSliderInput(session,inputId = "pu_new",
                              value = metadata$purity[metadata$sample == input$var])
        })

        output$original_fit <- shiny::renderPlot({
            orig_fit <- data.list[[input$var]]
            orig_purity <- metadata$purity[metadata$sample == input$var]

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

            orig_ploidy <- metadata$ploidy[metadata$sample == input$var]
            orig_purity <- metadata$purity[metadata$sample == input$var]

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

        output$fit_sunrise <- shiny::renderPlot({
            orig_ploidy <- metadata$ploidy[metadata$sample == input$var]
            orig_purity <- metadata$purity[metadata$sample == input$var]
            plotSunrise(data = data.list[[input$var]],
                        ploidy=orig_ploidy,purity=orig_purity)
        })

        output$fit_stat <- shiny::renderTable({
            fitTab <- calculateCINStats(data.list)
            fitTab <- as.data.frame(t(fitTab[rownames(fitTab) == input$var,]))
            return(fitTab)
        })

        output$new_stat <- shiny::renderTable({

            orig_ploidy <- metadata$ploidy[metadata$sample == input$var]
            orig_purity <- metadata$purity[metadata$sample == input$var]

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
            newTab <- as.data.frame(t(newTab[rownames(newTab) == input$var,]))
            return(newTab)
        })

        output$info <- shiny::renderTable({
            orig_ploidy <- metadata$ploidy[metadata$sample == input$var]
            orig_purity <- metadata$purity[metadata$sample == input$var]

            shiny::req(input$sunrise_click)
            #browser()
            shiny::nearPoints(calculateSunrise(data = data.list[[input$var]],
                                               ploidy=orig_ploidy,
                                               purity=orig_purity),
                       input$sunrise_click)[1,]
            #
            # x <- input$sunrise_click$x
            # y <- input$sunrise_click$y
            # cat("[", x, ", ", y, "]", sep = "")
        })
        output$qc_table <- DT::renderDataTable({
            DT::datatable(qcData$data,rownames = FALSE,extensions = 'Scroller',
                      options = list(pageLength = 5,
                                     ordering=F,
                                     scrollY = 200,
                                     scroller = TRUE,
                                     searchHighlight = T
                                     ))
        })

        output$saveQC <- shiny::downloadHandler(
            filename = function() {
                paste0(Sys.Date(),"_QC_table",".tsv")
            },
            content = function(file) {
                write.table(qcData$data,file,quote = F,sep = "\t",append = F,
                            row.names = F,col.names = T)
            }
        )

    }

    shiny::shinyApp(ui, server,options = list(launch.browser = TRUE))
}
