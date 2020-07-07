library(shiny)
library(shinydashboard)

ui <- dashboardPage(skin="yellow",   
  dashboardHeader(title = "The effect of missing data on design efficiency in cross-sectional multi-period parallel-arm cluster randomized trials.",disable = TRUE),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    titlePanel("The effect of missing data on design efficiency in cross-sectional multi-period parallel-arm cluster randomized trials"),
       fluidRow(
      column(width=4,
       
        box(
          title = "Sample size per cluster-period", width = NULL, solidHeader = FALSE, status = "primary",
          sliderInput("m", "Number of subjects per day (m)",min = 1, max = 100,value = c(1,10))
      ),
      
      tabBox(
        title = "Design specification",width=NULL,  
        id = "tabset1", 
        tabPanel("Design 1", 
                 checkboxGroupInput("checkDays1", 
                                    ("Select the days at which measurements are to be taken"), 
                                    choices = list("Monday" = 1, 
                                                   "Tuesday" = 2, 
                                                   "Wednesday" = 3,
                                                   "Thursday"=4,
                                                   "Friday"=5,
                                                   "Saturday"=6,
                                                   "Sunday"=7),
                                    selected = c(1,2,3,4,5,6,7),inline=TRUE),
                 numericInput("R1", "Number of weeks", 4, min = 1, max = 10),
                 numericInput("k1", "Number of clusters per condition", 10, min = 1, max = 100)
        ),
        tabPanel("Design 2", 
                 checkboxGroupInput("checkDays2", 
                                    ("Select the days at which measurements are to be taken"), 
                                    choices = list("Monday" = 1, 
                                                   "Tuesday" = 2, 
                                                   "Wednesday" = 3,
                                                   "Thursday"=4,
                                                   "Friday"=5,
                                                   "Saturday"=6,
                                                   "Sunday"=7),
                                    selected = c(1,2,3,4,5),inline=TRUE),
                 numericInput("R2", "Number of weeks", 4, min = 1, max = 10),
                 numericInput("k2", "Number of clusters per condition", 10, min = 1, max = 100)
        ),
        tabPanel("Design 3", 
                 checkboxGroupInput("checkDays3", 
                                    ("Select the days at which measurements are to be taken"), 
                                    choices = list("Monday" = 1, 
                                                   "Tuesday" = 2, 
                                                   "Wednesday" = 3,
                                                   "Thursday"=4,
                                                   "Friday"=5,
                                                   "Saturday"=6,
                                                   "Sunday"=7),
                                    selected = c(1,2,4,5),inline=TRUE),
                 numericInput("R3", "Number of weeks", 4, min = 1, max = 10),
                 numericInput("k3", "Number of clusters per condition", 10, min = 1, max = 100)
        ),
        tabPanel("Design 4", 
                 checkboxGroupInput("checkDays4", 
                                    ("Select the days at which measurements are to be taken"), 
                                    choices = list("Monday" = 1, 
                                                   "Tuesday" = 2, 
                                                   "Wednesday" = 3,
                                                   "Thursday"=4,
                                                   "Friday"=5,
                                                   "Saturday"=6,
                                                   "Sunday"=7),
                                    selected = c(1,2,3,4),inline=TRUE),
                 numericInput("R4", "Number of weeks", 4, min = 1, max = 10),
                 numericInput("k4", "Number of clusters per condition", 10, min = 1, max = 100)
        ),
        tabPanel("Design 5", 
                 checkboxGroupInput("checkDays5", 
                                    ("Select the days at which measurements are to be taken"), 
                                    choices = list("Monday" = 1, 
                                                   "Tuesday" = 2, 
                                                   "Wednesday" = 3,
                                                   "Thursday"=4,
                                                   "Friday"=5,
                                                   "Saturday"=6,
                                                   "Sunday"=7),
                                    selected = c(1,2,4),inline=TRUE),
                 numericInput("R5", "Number of weeks", 4, min = 1, max = 10),
                 numericInput("k5", "Number of clusters per condition", 10, min = 1, max = 100)
       ) )
      ),
      
      column(width=4,
             box(
               title = "Correlation parameters", width = NULL, solidHeader = FALSE, status = "primary",
               numericInput("ICC", "Intraclass correlation coefficient (rho)", 0.025, min = 0.001, max = 0.99),
               numericInput("decay", "Decay in variance per day (1-r)", 0.05, min = 0.001, max = 0.99)
             ),
             
      box(
        title = "Effect size and type I error rate", width = NULL, solidHeader = FALSE, status = "primary",
        numericInput("delta", "Standardized effect size (Cohen's d)", 0.2, min = 0, max = 0.8),
        numericInput("alpha", "Type I error rate (alpha)", 0.05, min = 0, max = 0.8),
        selectInput("test", label = "Type of test", 
                    choices = list("One-sided" = 1, "Two-sided" = 2), 
                    selected = 1))
      ),
      
      tabBox(
        title = "Attrition function",width=4,  
        # The id lets us use input$tabset1 on the server to find the current tab
        id = "tabset1", height = 500,
        tabPanel("Input parameters", 
                 
                 box(width=6,title="Control condition",
                 numericInput("omega0", "Parameter omega of Weibull survival function", 0.2, min = 0, max = 1),
                 numericInput("gamma0", "Parameter gamma of Weibull survival function", 2, min = 0, max = 100)),
           
                 box(width=6,title="Intervention condition",
                 numericInput("omega1", "Parameter omega of Weibull survival function", 0.2, min = 0, max = 1),
                 numericInput("gamma1", "Parameter gamma of Weibull survival function", 2, min = 0, max = 100)),
                 br(),br(),br(),
                 box(width=6,title="Scaling of time variable",
                     numericInput("maxweeks", "Maximum duration of the trial in weeks", 4, min = 0, max = 365))

                 ),
        tabPanel("Surival probability graph", 
                 fluidRow(
                   column(12,
                          plotOutput("survivalplot"))
                 )
                 ),
        tabPanel("Hazard probability graph", 
                 fluidRow(
                   column(12,
                          plotOutput("hazardplot"))
                 )
        )
      )
    ),

    fluidRow(  
   box( solidheader=FALSE,status="primary", width=4 , height = 550,
        br(),br(), 
        #img(src = "uu.png", height = 200, width = 200, align="center"),
        br(),br(),br(), br(), br(), 
      box(title = "Compare designs", width=12, solidHeader = TRUE, status = "warning",background = "yellow",
        submitButton("Submit")
      )),
      
      tabBox(
        title = "Output",width=8,  height = 550,
        # The id lets us use input$tabset1 on the server to find the current tab
        id = "tabset1", 
        tabPanel("Graphs", 
                 plotOutput("Resultsplot",width = "100%", height = "500px")
        ),
       
        tabPanel("Table", 
                 dataTableOutput(outputId = "ResultsTable")
        )
      )

    )
  )
)
