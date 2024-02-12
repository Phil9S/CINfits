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
        shiny::column(width = 6,
            shiny::fluidRow(
                shiny::column(width = 4,
                    shiny::selectInput("var", "Sample", choices = names(data.list)),
                    shiny::fluidRow(
                        shiny::column(6,shiny::checkboxInput("as",label = "allele-specific",value = AS)),
                        shiny::column(6,shiny::checkboxInput("round_values",label = "round segments",value = FALSE))
                        ),
                    shiny::conditionalPanel(condition = "input.as == true",
                                            shiny::wellPanel(
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
                                            )),
                    shiny::sliderInput("pl_new","ploidy",min = 1,max = 8,step = 0.01,value = 2),
                    shiny::sliderInput("pu_new","purity",min = 0.2,max = 1,step = 0.01,value = 0.7),
                    shiny::fileInput("qc_file",label = "QC file",multiple = F,
                                 accept = c(".tsv",".csv"),buttonLabel = "upload",
                                 placeholder = "upload existing qc file..."),
                    shiny::actionButton("accept_fit",label = "accept current"),
                    shiny::actionButton("refit_fit",label = "accept refit"),
                    shiny::actionButton("reject_fit",label = "reject"),
                    shiny::textAreaInput("notes_fit",label = "notes",
                                         placeholder = "Add fit notes here...")
                ),
                shiny::column(width = 8,
                    shiny::plotOutput("fit_sunrise",click = "sunrise_click"),
                    shiny::tableOutput(outputId = "info")
                ),
            ),
            shiny::fluidRow(
                shiny::column(width = 12,
                    DT::dataTableOutput(outputId = "qc_table")
                )
            ),
        ),
        shiny::column(width = 6,
            shiny::fluidRow(
                shiny::column(width = 12,
                    shiny::fluidRow(
                        shiny::column(width = 9,
                            shiny::plotOutput("original_fit")),
                        shiny::column(width = 2,
                            shiny::tableOutput("fit_stat"))
                    ),
                    shiny::fluidRow(
                    shiny::column(width = 9,
                        shiny::plotOutput("new_fit")),
                    shiny::column(width = 2,
                        shiny::tableOutput("new_stat"))
                    )
                )
            )
        )
    )

    server <- function(input, output, session) {

        # get allele-specific status and modulate checkbox
        shiny::observe(if(input$as == TRUE | input$as == FALSE){
            shinyjs::disable(id = "as")
        })

        # setting sample index position
        samplePos = shiny::reactiveVal(1)

        # observe QC table
        qcData <- shiny::reactiveValues(data=NULL)
        shiny::observe({
            qcData$data <- qctable
        })
        # action on accepting fit
        shiny::observeEvent(input$accept_fit,{
            samplePos(samplePos()+1)
            selectedSample <- names(data.list)[samplePos()]

            df <- qcData$data
            df$use[df$sample == input$var] <- TRUE
            qcData$data <- df

            shiny::updateSelectInput(inputId = "var",selected = selectedSample)

        })

        # Update sliders to provided fit
        shiny::observe({
            shiny::updateSliderInput(inputId = "pl_new",
                              value = metadata$ploidy[metadata$sample == input$var])
            shiny::updateSliderInput(inputId = "pu_new",
                              value = metadata$purity[metadata$sample == input$var])
        })

        output$original_fit <- shiny::renderPlot({
            orig_fit <- data.list[[input$var]]
            orig_purity <- metadata$purity[metadata$sample == input$var]
            if(input$round_values){
                orig_fit$segVal <- round(orig_fit$segVal)
                if(AS){
                    orig_fit$nAraw <- round(orig_fit$nAraw)
                    orig_fit$nBraw <- round(orig_fit$nBraw)
                }
            }

            if(AS){
                cols <- c(input$nA,input$nB)
            } else {
                cols <- "red"
            }

            plotProfile(orig_fit,sample = input$var,
                        cn.max = 15,purity = orig_purity,
                        alleleSpecific = AS,cols = cols)
        })
        output$new_fit <- shiny::renderPlot({

            orig_ploidy <- metadata$ploidy[metadata$sample == input$var]
            orig_purity <- metadata$purity[metadata$sample == input$var]

            new_fit <- rescaleFit(data = data.list[[input$var]],
                                  old_ploidy = orig_ploidy,
                                  old_purity = orig_purity,
                                  new_ploidy = input$pl_new,
                                  new_purity = input$pu_new)

            if(input$round_values){
                new_fit$segVal <- round(new_fit$segVal)
                if(AS){
                    # new_fit$nAraw <- round(new_fit$nAraw)
                    # new_fit$nBraw <- round(new_fit$nBraw)
                }
            }

            if(AS){
                cols <- c(input$nA,input$nB)
            } else {
                cols <- "red"
            }

            plotProfile(new_fit,sample = input$var,
                        cn.max = 15,purity = input$pu_new,
                        alleleSpecific = AS,cols = cols)
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
                                 new_purity = input$pu_new)
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
                                     #lengthMenu = c(5, 20, 50, 100),
                                     searchHighlight = T#,
                                     #deferRender = TRUE,
                                     #scrollY = 200,
                                     #scroller = TRUE)
                                     ))
        })

    }

    shiny::shinyApp(ui, server,options = list(launch.browser = TRUE))
}
