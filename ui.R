library(shiny)
library(fpow)
shinyUI(fluidPage(
  titlePanel(h1("Twin power calculator")),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("nmz",
               label = h3("Number of MZ pairs"),
               value = 20, min = 0.0, step = 5),
      numericInput("ndz",
               label = h3("Number of DZ pairs"),
               value = 20, min = 0.0, step = 5),
      numericInput("Va",
               label = h3("Additive genetic variance"),
               value = 0.6, min = 0.0, step = 0.1),
      numericInput("Vc",
               label = h3("Common environmental variance"),
               value = 0.2, min = 0.0, step = 0.1),
      numericInput("Ve",
               label = h3("Environmental variance"),
               value = 0.2, min = 0.0, step = 0.1),
      numericInput("alpha",
               label = h3("Type-1 error rate"),
               value = 0.05, min = 0.0, max = 0.499999, step = 0.05),
      numericInput("power",
               label = h3("Required power"),
               value = 0.8, min = 0.0, max = 1.0, step = 0.05),
      actionButton("goButton", "Calculate!")
      ),
    mainPanel(
         h1("This Shiny application provides automated power analysis for 
             the detection of additive genetic (A) and common environmental 
             (C) variance components of a quantitative trait in the 
             classical twin design. "),
         br(),
         h2("Adjust the parameters on the left and press calculate to view the results"), 
         br(),
      htmlOutput("text1"),
      h2("Results"),
      tableOutput("view1"),
      fluidRow(column(width = 5,
                   strong("ACE vs. CE model"),
                   tableOutput("view2")),
            column(width = 5,
                   strong("ACE vs. AE model"),
                   tableOutput('view3'))
       ),
       br(),
         br(),
         br(),
         br(), 
       strong("Definitions"),
       p("NCP: Noncentrality parameter, 
          ML:  Maximum likelihood, 
          E:   Environmental variance component,
          MZ:  Monozygotic,
          DZ:  Dizygotic"),
       strong("If you use this site, please reference the following:"),
         a("Visscher P.M. (2004). 'Power of the classical twin design revisited'. 
           Twin Research 7, 505-512 and Visscher P.M., Gordon S., Neale M.C. (2008) 
           and 'Power of the classical twin design revisited: II Detection of Common 
           Environmental Variance'. Twin Research and Human Genetics 11, 48-54."),
       p(" "),
       p(" "), 
       p("Please send suggestions or comments to Peter Visscher or if you are having
            difficulties with the application please contact Luke Lloyd-Jones at
            l.lloydjones@uq.edu.au. This application was prepared by Matthew Robinson,
            Luke Lloyd-Jones and Peter Visscher.")
  ))
  ))