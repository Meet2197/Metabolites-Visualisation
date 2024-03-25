install.packages("shiny")
install.packages("shinydashboard")
library(shiny)
library(shinydashboard)

# Define UI for application
ui <- fluidPage(
  tags$style(HTML('.btn-tab { 
                    display: inline-block; 
                    margin: 5px; 
                    background-color: #007bff; 
                    color: #fff; 
                    border: none; 
                    padding: 10px 20px; 
                    border-radius: 5px; 
                    cursor: pointer; 
                    font-size: 16px; 
                    text-align: center; 
                    text-decoration: none;
                    }
                  .btn-tab.active { 
                    background-color: #0056b3; 
                  }
                ')),
  navbarPage(
    title = "Dashboard with Tabs",
    tabPanel("Tab 1",
             fluidRow(
               column(6, textInput("text1", "Text Box 1")),
               column(6, textInput("text2", "Text Box 2")),
               column(6, textInput("text3", "Text Box 3")),
               column(6, textInput("text4", "Text Box 4"))
             )
    ),
    tabPanel("Tab 2",
             fluidRow(
               column(6, textInput("text5", "Text Box 5")),
               column(6, textInput("text6", "Text Box 6")),
               column(6, textInput("text7", "Text Box 7")),
               column(6, textInput("text8", "Text Box 8"))
             )
    ),
    tabPanel("Tab 3",
             fluidRow(
               column(6, textInput("text9", "Text Box 9")),
               column(6, textInput("text10", "Text Box 10")),
               column(6, textInput("text11", "Text Box 11")),
               column(6, textInput("text12", "Text Box 12"))
             )
    ),
    tabPanel("Tab 4",
             fluidRow(
               column(6, textInput("text13", "Text Box 13")),
               column(6, textInput("text14", "Text Box 14")),
               column(6, textInput("text15", "Text Box 15")),
               column(6, textInput("text16", "Text Box 16"))
             )
    ),
    position = "static-top"
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("user", "Select User", c("User 1", "User 2", "User 3"))
    ),
    mainPanel(
      
    )
  )
)

# Define server logic
server <- function(input, output) {
  
}

# Run the application
shinyApp(ui = ui, server = server)
  

# Run the application
