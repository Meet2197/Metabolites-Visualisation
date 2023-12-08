library(shiny)
library(plotly)
library(DescTools)
library(dplyr)
shinyServer(function(input, output) {
    data <- reactive({  
    dist <- switch(input$dist,
                   norm = rnorm,
                   unif = runif,
                   lnorm = rlnorm,
                   exp = rexp,
                   rnorm)
    dist(input$n)
  })
  output$plot <- renderPlot({
    dist <- input$dist
    n <- input$n

    hist(data(),
         main=paste('r', dist, '(', n, ')', sep=''))
    
    })
  
  output$plot1 <- renderPlotly({
    f.freq <- Freq(f$SAMPLETYPE)
    plot_ly(f.freq, x = f.freq$freq, y = f.freq$level, color = f.freq$level, type = "bar") %>% layout(title = "Samples Variation Clustering")
  })
  output$summary <- renderPrint({
    summary(data())
  })
  output$table <- renderTable({
    data.frame(f.freq[1:12, ])
  })
})
