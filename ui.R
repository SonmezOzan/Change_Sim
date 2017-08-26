library(shiny)
library(shinydashboard)
library(fda)
library(plotly)


ui <- dashboardPage(
   dashboardHeader(title="Simulation"
   ),
   ## Sidebar content
   dashboardSidebar(
      sidebarMenu(
         menuItem("Home", tabName = "home", icon = icon("home")),
         menuItem("Analysis", tabName = "widgets", icon = icon("th")),
         menuItem("Contact", tabName = "contact", icon = icon("dashboard"))
      )
   ),
   dashboardBody(
      tabItems(
         tabItem(tabName = "home",
                 h1("Change Point Analysis in Functional Data"),
                 h4("This interactive application enables interested audiance to explore
                    change point analysis in the mean function of functional data utilizing
                    fPCA-based methodology and the developed fully functional-methodology."),
                 h4("This simulation app is a supplement to the below paper "),
                 wellPanel(
                    helpText(   a("Detecting and dating structural breaks in functional data without dimension reduction",  
                                  href="https://arxiv.org/pdf/1511.04020v2.pdf")
                    )),
                 h4("The ideas that are presented in this paper are also applied to Australian
                        temperature data that are recoreded from 8 different stations. The change point
                    analysis of Australian temperature can be found in the below link."),
                 wellPanel(
                    helpText(a("Australian Temperature Data Analysis Shiny App",  
                               href="https://changepointappozan.shinyapps.io/AustraliaTemperatureData/")
                    )),
                    h4("We strongly encourage to read the manuscript before exploring the simulatio App,
                       as well as the Australian Temperature Data Analysis App")
                 ),
         # Second tab content
         tabItem(tabName = "widgets",
                 sidebarLayout(
                    sidebarPanel(
                    selectizeInput('N','Sample Size (n)', choices = c(50, 100, 200, 500, 1000)),
                    radioButtons("DataK", "Data Type",
                                 c("i.i.d." = "IID",
                                   "FAR(1)" = "FAR")),
                    radioButtons("Sigma", "Setting for sigma",
                                 c("Fast Setting" = "fast",
                                   "Slow Setting" = "slow",
                                   "Custom Setting" = "custom")),
                    sliderInput("k", 
                                "Change Function (function of k)", 
                                value = 1,
                                min = 1, 
                                max = 21),
                    sliderInput("theta", 
                                "Change location (theta)", 
                                value = 0.5,
                                min = 0, 
                                max = 1),
                    numericInput("SNR", "Signal-to-Noise Ratio (SNR)", 1),
                    radioButtons("TVE", "TVE for fPCA-based method", c(0.85, 0.9, 0.95)),
                    submitButton("Update View")),
                 mainPanel(
                    plotOutput("plotd"),
                    tableOutput("ff"),
                    plotlyOutput("plot")
                 ))
                 
         ),
         tabItem(tabName = "contact", 
                 h3("Ozan Sonmez"),
                 h5("Department of Statistics"),
                 h5("Mathematical Sciences Building 1214"),
                 h5("University of California, Davis"),
                 h5("One Shields Avenue"),
                 h5("Davis, CA 95616"),
                 h4("osonmez@ucdavis.edu"))
      )
   )
)

