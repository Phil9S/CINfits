#' fitInteractive
#'
#' Launch an interactive shiny application  to perform copy number fitting and
#' quality control
#'
#' @param data data.frame or list of segmented copy number profiles
#'
#' @return interactive shiny application
#' @export
#'
fitInteractive <- function(data=NULL){
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop(
            "Package \"shiny\" must be installed to use interactive fitting",
            call. = FALSE
        )
    }
    if(is.null(data)){
        stop("no data")
    }
    data.list <- split(data,f = data$sample)
    ui <- shiny::fluidPage(
        title = "CNfits",
        shiny::titlePanel(title = "CNfits"),
        shiny::fluidRow(
            shiny::column(width = 2,
                shiny::selectInput("var", "Variable", choices = names(data.list)),
                shiny::sliderInput("pl_new","ploidy",min = 1,max = 8,step = 0.1,value = 2),
                shiny::sliderInput("pu_new","purity",min = 0.2,max = 1,step = 0.01,value = 0.7),
                shiny::checkboxInput("round_values",label = "round segments",value = FALSE)),
            shiny::column(width = 4,
                shiny::plotOutput("fit_sunrise",click = "sunrise_click"),
                shiny::tableOutput(outputId = "info")),
            shiny::column(width = 4,
                shiny::plotOutput("original_fit"),
                shiny::plotOutput("new_fit"))
            ),
        shiny::fluidRow(
            shiny::column(width = 4,shiny::tableOutput("fit_stat")),
            shiny::column(width = 4,shiny::tableOutput("new_stat"))
        )
    )

    server <- function(input, output, session) {
        output$original_fit <- shiny::renderPlot({
            orig_fit <- data.list[[input$var]]

            if(input$round_values){
                orig_fit$segVal <- round(orig_fit$segVal)
            }

            plotProfile(orig_fit,sample = input$var,cn.max = 15)
        })
        output$new_fit <- shiny::renderPlot({
            orig_ploidy <- calculatePloidy(data.list[[input$var]])
            orig_purity <- 0.7 # TEMP
            new_fit <- rescaleFit(data = data.list[[input$var]],
                                  ploidy = input$pl_new,purity = input$pu_new)

            if(input$round_values){
                new_fit$segVal <- round(new_fit$segVal)
            }

            plotProfile(new_fit,sample = input$var,cn.max = 15)
        })
        output$fit_sunrise <- shiny::renderPlot({
            plotSunrise(data = data.list[[input$var]])
        })

        output$fit_stat <- shiny::renderTable({
            fitTab <- calculateCINStats(data.list)
            fitTab <- as.data.frame(t(fitTab[rownames(fitTab) == input$var,]))
            return(fitTab)
        })

        output$new_stat <- shiny::renderTable({
            newTab <- rescaleFit(data = data.list[[input$var]],
                                  ploidy = input$pl_new,purity = input$pu_new)
            newTab <- calculateCINStats(newTab)
            newTab <- as.data.frame(t(newTab[rownames(newTab) == input$var,]))
            return(newTab)
        })

        output$info <- shiny::renderTable({
            shiny::req(input$sunrise_click)
            #browser()
            shiny::nearPoints(calculateSunrise(data = data.list[[input$var]]),
                       input$sunrise_click)[1,]
            #
            # x <- input$sunrise_click$x
            # y <- input$sunrise_click$y
            # cat("[", x, ", ", y, "]", sep = "")
        })
    }

    shiny::shinyApp(ui, server)
}
