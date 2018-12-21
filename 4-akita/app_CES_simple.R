# CES with production model
pkg.name <- "shiny"
if(!require(pkg.name, character.only=TRUE)){
  install.packages(pkg.name)
}

library(shiny)

# Define UI for application that draws a histogram
ui <-  fluidPage(
  titlePanel("産卵親魚量（取り残し）一定方策 Ver.12/20"),
  
  fluidRow(
    column(2,
           wellPanel(
             sliderInput("r","増殖率:",min = 0,max = 2,value = 1,step = 0.1),
             sliderInput("escape","取り残し量:",min = 0,max = 50,value = 30,step = 10),
             sliderInput("sigma","環境変動:",min = 0, max = 1,value = 0,step = 0.1),
             sliderInput("num","シミュレーション回数:",min = 1,max = 10,value = 1)
           )),
    
    column(10,plotOutput("plot_N"))
  )
)

server <- function(input, output) {
  
  output$plot_N <- renderPlot({
    r <- input$r
    K <- 100
    num <- input$num
    N <- matrix(0,100,num)
    catch <- matrix(0,100,num)
    sigma <- input$sigma
    escape <- input$escape
    N[1,] <- K
    for(i in 2:100){
      for (j in 1:num){
        catch[i,j] <- ifelse((N[i-1,j] - escape) > 0, N[i-1,j] - escape, 0)
      }
      N[i,] <- (N[i-1,] + r*N[i-1,]*(1-N[i-1,]/K) - catch[i,])*exp(rnorm(num,-0.5*sigma^2,sigma))
      N[i,N[i,]<0] <- 0
      
    }
    par(mfrow=c(2,1))
    par(oma = c(2, 2, 0, 0))
    par(mar = c(3, 5, 1, 1))
    matplot(N,ylab="Biomass", ylim=c(0,100),type="l",lty=1,cex.lab=2,cex.axis=2,lwd=2)
    matplot(catch,ylab="Catch", ylim=c(0,max(catch)),type="l",lty=1,cex.lab=2,cex.axis=2,lwd=2)
    mtext("Year",side = 1,    # 1:下側、2:左側、3:上側、4:右側
          line = 1,               # マージン領域の何行目か    
          outer = TRUE,cex = 2
          )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

