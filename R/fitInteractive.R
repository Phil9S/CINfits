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
    ui <- fluidPage(
            column(width = 4,
                selectInput("var", "Variable", choices = names(data.list)),
                sliderInput("pl_new","ploidy",min = 1,max = 8,step = 0.1,value = 2),
                sliderInput("pu_new","purity",min = 0.2,max = 1,step = 0.01,value = 0.7)),
            column(width = 4,
                plotOutput("fit_sunrise")),
            column(width = 4,
                plotOutput("original_fit"),
                plotOutput("new_fit")
            )
    )

    server <- function(input, output, session) {
        output$original_fit <- renderPlot({
            plotprofile(data.list[[input$var]],sample = input$var,cn.max = 15)
        })
        output$new_fit <- renderPlot({
            orig_ploidy <- calculatePloidy(data.list[[input$var]])
            orig_purity <- 0.7 # TEMP
            new_fit <- rescaleFit(data = data.list[[input$var]],ploidy = input$pl_new,purity = input$pu_new)
            plotprofile(new_fit,sample = input$var,cn.max = 15)
        })
        output$fit_sunrise <- renderPlot({
            plotSunrise(data = data.list[[input$var]])
        })
    }

    shinyApp(ui, server)
}
