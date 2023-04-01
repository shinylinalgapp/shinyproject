library(shiny)
library(ggplot2)
library(fastmatrix)
library(microbenchmark)
source("matrixalgs.R")

# Define UI for inputs

inputPanel <-  sidebarPanel(
  sliderInput("n", "Matrix size n", min = 2, max = 200, value = 100),
  
  br(),
  
  # SOR parameter w
  numericInput("w",
               "SOR parameter w",
               value = 1.5,
               min = 0,
               max = 2,
               step = 0.01),
  numericInput("sparsity",
               "Sparsity",
               value = 0,
               min = 0,
               max = 1,
               step = 0.01),
  numericInput("tol", "Error Tolerance",
               value = 10^(-6),
               min = 10^(-12),
               max = 10^(-2),
               step = 10^(-12)),
  numericInput("dist", "Initial Approximation Error",
               value = 100,
               min = 0,
               max = 5000,
               step = 1),
  sliderInput("num", "Number of equations per coefficient matrix", 
              value = 1,
              min = 1, 
              max = 10),
  actionButton("run", "Run Comparison"))
ui <- fluidPage(
  
  # App title ----
  titlePanel("Intro to Numerical Linear Algebra"),

  navbarPage("",
                  tabPanel("About", 
                           "This app is a tool for students to explore various numerical algorithms for solving systems of linear equations.",
                           "This problem comes up extensively in applications such as balancing chemical equations in chemistry, as well as numerous other applications in physics and engineering.",
                           br(),
                           br(),
                           "The typical algorithm students learn in an introductory linear algebra course for solving systems of linear equations is Gaussian Elimination.",
                           "This algorithm consists of repeatedly reducing rows using basic row operations to find an equivalent upper triangular matrix.",
                           "As intuitive as the algorithm is, it isn't used in practice; the LU factorization is a modified version of Gaussian Elimination",
                           "which factors a matrix by dynamically storing the factors used in the row operations of Gaussian Elimination.",
                           "This factorization allows many problems with the same coefficient matrix to be solved in quadratic time, after an initial cubic time factorization.",
                           br(),
                           br(),
                           "Alternatively, iterative methods exist like Gauss-Seidel which continuously improve an initial approximation,",
                          "rather than finding an exact answer. These iterative methods work very well when the coefficient matrix is sparse; that is, it has few non-zero elements.",
                          "SOR (Successive Over Relaxation) is a modification of Gauss Seidel that, on each iteration,",
                          "computes a weighted average of the current approximation and the next approximation from the Gauss-Seidel algorithm. Depending on the value of the weight w,",
                          "this can result in faster convergence rates than Gauss-Seidel. When w = 1, the methods are identical.",
                          "However, unlike LU factorization and Gaussian Elimination, these methods do not work for all invertible coefficient matrices;",
                          "they may fail to converge if the matrix is not strictly diagonally dominant. A matrix is strictly diagonally dominant if, for each row,",
                          "the magnitude of the element that is on the diagonal is greater than the sum of the magnitudes of all other elements.",
                          br(),
                          br(),
                          "Therefore, there is an important algorithmic choice to be made regarding which algorithm to choose when in one of the numerous applications they are used in.",
                          "This app is designed to assist in making that choice in a fun and interactive manner. The Algorithm Comparison tab contains interactive graphs to play around with ",
                          "which are dynamic and their parameters can be modified in the sidebar panel."),
                  tabPanel("Algorithm Comparison", sidebarLayout(inputPanel, mainPanel("Press Run Comparison to generate graphs!",
                                                                                       plotOutput("runtimeComp"),
                                                                                       plotOutput("seidelPlot"),
                                                                                       plotOutput("sorPlot")))),
             tabPanel("Design Process",
                      "This app was designed while consulting T.Moran's system design book and the essence of software book by Daniel Jackson. The latter provided very useful information ",
                      "about concepts and how to design them well, which is reflected in the app. The concepts of the app are very clear and defined, and are generally a one-to-one ",
                      "correspondance with features on the app. For example, the matrix size concept is very clearly defined and is given a numeric input ",
                      "that can be directly modified by the user in an unambiguous way. The other concepts of the app, like the sparsity concept, behave similarly.",
                      "As such, the use of concepts greatly informed the design process, and the app was designed in a manner that makes them as clear as possible ",
                      "for a good user experience.", br(), br(), 
                      "T. Moran's book also very greatly influenced the design of the app. I was careful to make sure the app was very easy to learn, easy to use once learned, ",
                      "and limit the number of possible errors the user can make."),
             tabPanel("Acknowledgments and References", 
                      h1("Acknowledgements"), "I would like to thank my CS professor for providing continuous guidance and feeback on my project idea.",
                      h1("References"), 
                      h2("Numerical Algorithms References"), 
                      "This app was heavily inspired by the material from Chapter 2 of Timothy Sauer's 2013 Numerical Analysis textbook.", 
                      "Additionally, the Wikipedia pages for each of the algorithms were used when designing the app and implementing the algorithms.",
                      h2("R references"),
                      "The R Shiny Tutorial (https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/) was extensively used throughout the creation of this app,",
                      "as was the Engineering Production Grade Shiny apps book.",
                      "Additionally, the documentation for the packages ggplot2, fastmatrix, and microbenchmark were frequently consulted.",
                      "The tabsets R shiny example was used as a base for this project.",
                      h2("Design references"),
                      "The designs of R shiny apps in the Shiny user showcase (https://shiny.rstudio.com/gallery/) were consulted extensively.",
                      "In addition, T.Moran's system design book and Daniel Jackson's the essence of software book were also used.")))
                  
      
      
    

server <- function(input, output) {
  
  # Runs and times gausselim using provided parameters
  gaussTime <- eventReactive(input$run, {
    gaussTimes <- 1:10
    for (i in 1:10) {
      matrix <- generateSddMatrix(input$sparsity, input$n)
      b  <- sample(vals, input$n, replace = TRUE)
      gaussTimes[i] <- microbenchmark::microbenchmark(gaussian_elimination(matrix, b, input$n), times = 1)$time
    }
    gaussTimes
  })
  
  sorTime <- eventReactive(input$run, {
    sorTimes <- 1:10
    for (i in 1:10) {
      matrix <- generateSddMatrix(input$sparsity, input$n)
      vector  <- sample(vals, input$n, replace = TRUE)
      vector2 <- generateVector(gaussian_elimination(matrix, vector, input$n), input$dist, input$n)
      sorTimes[i] <- microbenchmark::microbenchmark(SOR(input$w, matrix, vector, vector2, input$tol, 1000), times = 1)$time
    }
    sorTimes
  })
  
  seidelTime <- eventReactive(input$run, {
    seidelTimes <- 1:10
    errors <- 1:10
    for (i in 1:10) {
      matrix <- generateSddMatrix(input$sparsity, input$n)
      vector  <- sample(vals, input$n, replace = TRUE)
      vector2 <- generateVector(gaussian_elimination(matrix, vector, input$n), input$dist, input$n)
      seidelTimes[i] <- microbenchmark::microbenchmark(SOR(1, matrix, vector, vector2, input$tol, 1000), times = 1)$time
    }
    seidelTimes
  })
  
  luTime <- eventReactive(input$run, {
    luTimes <- 1:10
    for (i in 1:10) {
      matrix <- generateSddMatrix(input$sparsity, input$n)
      vector  <- sample(vals, input$n, replace = TRUE)
      luTimes[i] <- microbenchmark::microbenchmark(A <- lu(matrix), times = 1)$time
      sum <- 0
      for (j in 1:input$num)
      {
        sum <- sum + microbenchmark::microbenchmark(solve.lu(A, vector), times = 1)$time
      }
      luTimes[i] <- (sum*10 + luTimes[i])/(input$num)
    }
    luTimes
  })
  
  # Display the average time in the UI using renderText
  output$gaussAvg <- renderText({
    if (input$run == 0) {
      # Default 0
      0
    } else {
      # Otherwise average gauss time
      gaussTime()
    }
  })
  output$sorAvg <- renderText({
    if (input$run == 0) {
      # Default 0
      0
    } else {
      # Otherwise average SOR time
      sorTime()
    }
  })
  output$seidelAvg <- renderText({
    if (input$run == 0) {
      # Default 0
      0
    } else {
      # Otherwise average seidel time
      seidelTime()
    }
  })
  
  
  data <- reactive({
    data.frame(
      Algorithm = c("Gaussian Elimination", "LU", "Gauss-Seidel", "SOR"),
      Runtime = c(mean(gaussTime()), mean(luTime()), mean(seidelTime()), mean(sorTime()))
    )
  })
  
  output$runtimeComp <- renderPlot({
    ggplot(data(), aes(x = Algorithm, y = Runtime)) +
      geom_bar(stat = "identity", fill = c("red", "purple", "blue", "green"))
  })
seidelErr <- eventReactive(input$run, {
    matrix <- generateSddMatrix(input$sparsity, input$n)
    vector  <- sample(vals, input$n, replace = TRUE)
    vector2 <- generateVector(vector, input$dist, input$n)
    SOR(1, matrix, vector, vector2, input$tol, 1000)
})

sorErr <- eventReactive(input$run, {
  matrix <- generateSddMatrix(input$sparsity, input$n)
  vector  <- sample(vals, input$n, replace = TRUE)
  vector2 <- generateVector(vector, input$dist, input$n)
  SOR(input$w, matrix, vector, vector2, input$tol, 1000)
})

seidelData <- reactive({
  data.frame(
    Error = seidelErr(),
    SeidelIterations = (1:length(seidelErr()))
    
  )
})

sorData <- reactive({
  errors <- sorErr()
  data.frame(
    SORIterations = (1:length(errors)),
    Error = errors
  )
})

output$seidelPlot <- renderPlot({
  ggplot(seidelData(), aes(x = SeidelIterations, y = Error)) +
    geom_line(color = "blue")
})

output$sorPlot <- renderPlot({
  ggplot(sorData(), aes(x = SORIterations, y = Error)) +
    geom_line(color = "green")
})


}


# Create Shiny app ----
shinyApp(ui, server)