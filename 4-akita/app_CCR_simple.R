# CCR with production model
pkg.name <- "shiny"
if(!require(pkg.name, character.only=TRUE)){
  install.packages(pkg.name)
}

library(shiny)

# Define UI for application that draws a histogram
ui <-  fluidPage(
  titlePanel("漁獲量一定方策 Ver.12/20"),
  
  fluidRow(
    column(2,
           wellPanel(
             sliderInput("r","増殖率:",min = 0,max = 2,value = 1,step = 0.1),
             sliderInput("catch","年漁獲量:",min = 0,max = 50,value = 20,step = 10),
             sliderInput("sigma","環境変動:",min = 0, max = 1,value = 0,step = 0.1),
             sliderInput("num","シミュレーション回数:",min = 1,max = 10,value = 1)
           )       
    ),
    
    column(10,
           plotOutput("plot_N", width = "100%")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot_N <- renderPlot({
    r <- input$r
    K <- 100
    num <- input$num
    N <- matrix(0,100,num)
    sigma <- input$sigma
    catch <- input$catch
    N[1,] <- K
    for(i in 2:100){
      N[i,] <- (N[i-1,] + r*N[i-1,] * (1-N[i-1,] / K) - catch) * exp(rnorm(num,-0.5*sigma^2,sigma))
      N[i,N[i,]<0] <- 0
    }
    par(mar = c(5, 5, 1, 1))
    matplot(N,xlab="Year",ylab="Biomass", ylim=c(0,K),type="l",lty=1,cex.lab=2,cex.axis=2,lwd=2)
  })
}

shinyApp(ui = ui, server = server)