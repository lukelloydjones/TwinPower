library(shiny)
library(fpow)
shinyUI(fluidPage(
  titlePanel("TwinPower calculator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("nmz",
               label = h3("Number of MZ pairs"),
               value = 20),
      numericInput("ndz",
               label = h3("Number of DZ pairs"),
               value = 20),
      numericInput("Va",
               label = h3("Additive genetic variance"),
               value = 0.6),
      numericInput("Vc",
               label = h3("Common environmental variance"),
               value = 0.2, min = 0.0, max = 1.0),
      numericInput("Ve",
               label = h3("Environmental variance"),
               value = 0.2, min = 0.0, max = 1.0),
      numericInput("alpha",
               label = h3("Type-1 error rate"),
               value = 0.05, min = 0.0, max = 0.499999),
      numericInput("power",
               label = h3("Required power"),
               value = 0.8, min = 0.0, max = 1.0),
      actionButton("goButton", "Calculate!")
      ),
    mainPanel(
     tabsetPanel(id = "tabset",
         h2("This Shiny app performs power calculations for variance components
             estimated under model --- using data from monozygotic and dizygotic 
             twins. Alter the parameters on the left panel and press calculate 
             to observe the results.")),
      htmlOutput("text1")
      )
  )
  ))