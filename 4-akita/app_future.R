# 漁業制御ルールの実践
# 推定は一切していない点に注意
# （推定誤差・バイアスは無視）
# 管理下では、資源量が管理水準を下回った場合にFを管理する

pkg.name <- "shiny"
if(!require(pkg.name, character.only=TRUE)){
  install.packages(pkg.name)
}

library(shiny)

ui <-  fluidPage(
  titlePanel("管理基準値に基づいた管理・将来予想 Ver.12/25"),
  
  fluidRow(
    column(2,helpText("漁獲係数：")),
    column(1,checkboxInput("Fmax", "Fmax", TRUE)),
    column(1,checkboxInput("F01", "F0.1", FALSE)),
    column(2,checkboxInput("FSPR20", "F20%SPR", FALSE)),
    column(2,checkboxInput("FSPR30", "F30%SPR", FALSE)),
    column(2,checkboxInput("FSPR40", "F40%SPR", FALSE)),
    #column(1,checkboxInput("Fmed", "Fmed", FALSE)),
    column(1,checkboxInput("Fmsy", "Fmsy", FALSE))
  ),
  
  fluidRow(
    column(2,helpText("資源水準：")),
    column(2,checkboxInput("SSB20", "20%SSB0", TRUE)),
    column(1,checkboxInput("R50", "R50%", FALSE)),
    column(1,checkboxInput("Bmsy", "Bmsy", FALSE))
  ),
  
  fluidRow(
    column(2,
           wellPanel(
             sliderInput("maturity_age","成熟齢:",min = 0,max = 10,value = 4,step = 1),
             sliderInput("catch_start_age","漁獲開始齢:",min = 0,max = 10,value = 0,step = 1),
             sliderInput("F_value","漁獲圧:",min = 0,max = 1,value = 0.35,step = 0.05),
             sliderInput("steepness","スティープネス:",min = 0.25,max = 1,value = 0.8,step = 0.05),
             sliderInput("sigmaR","加入変動:",min = 0, max = 1,value = 0,step = 0.1),
             sliderInput("sigma","努力量変動:",min = 0, max = 1,value = 0,step = 0.1)
           )),
    
    column(10,plotOutput("plot_N"))
  )
)

server <- function(input, output) {
  options(digits=4) # 出力を見やすくするため
  output$plot_N <- renderPlot({
    par(mfrow=c(2,4))
    par(oma = c(0, 0, 0, 0))
    par(mar = c(5, 5, 1.5, 1))
    par(cex.axis=2,cex.lab=2,cex.main=2,lwd = 2)
    F_RP <- B_RP <- NULL
    age <- 10 # 最大年齢クラス（この場合０〜９歳）
    maturity_age <- input$maturity_age # 成熟年齢
    M_at_age <- rep(0.3,age) # 自然死亡率は齢に依らず一定と仮定 
    Q_at_age <- c(rep(0,maturity_age),rep(1,age-maturity_age)) # 成熟年齢に達すると100%成熟すると仮定
    F_at_age <- c( rep(0,input$catch_start_age), rep(input$F_value, age - input$catch_start_age))
    age_vector <- seq(0,(age-1),1)
    #################################
    VB_a0 <- -0.1
    VB_b <- 3
    VB_K <- 0.5
    VB_winf <- 10
    
    VB <- function(a){
      VB_winf*(1 - exp(-VB_K*(a-VB_a0)) )^VB_b 
    }
    
    W_at_age <- sapply(seq(0,age-1,1),VB)
    
    BH_a <- 1
    BH_b <- 1
    
    BH <- function(SSB,BH_a,BH_b){
      R <- BH_a * SSB / (BH_b + SSB)
      return(R)
    }
    # no fishing
    burn_in <- 10
    SSB_init <- 5
    SSB <- c(SSB_init,rep(0,burn_in))
    N_at_age <- rep(0,age)
    for (y in 1:burn_in) {
      for (a in 1:age){
        if(a==1) {
          N_at_age[a] <- BH(SSB[y],BH_a,BH_b)  
        }else{
          N_at_age[a] <- N_at_age[a-1] * exp(-M_at_age[a-1])
        }
      }
      SSB[y+1] <- sum(N_at_age*Q_at_age*W_at_age)
    }
    S_at_age <- c(1,rep(0,age-1))
    for (a in 1:(age-1)){
      S_at_age[a+1] <- S_at_age[a] *exp(-M_at_age[a])
    }
    SPR0 <- sum(S_at_age * W_at_age * Q_at_age)
    #curve(expr =  BH_a * x / (BH_b + x),from = 0,to=8)
    #curve(expr = x/SPR0,add = T)
    B0 <- SSB_init2 <- SSB[y+1]
    R0 <- BH(SSB_init2,BH_a,BH_b)
    
    BH2 <- function(SSB,B0,R0,steepness){#reparameterization
      R <- 4*steepness*R0*SSB / (B0*(1-steepness) + SSB*(5*steepness-1))
      return(R)
    }
    
    #################################
    # fishing
    term <- 50
    sd <- input$sigmaR#approximatelly CV of Recruitiment
    sd2 <- input$sigma
    steepness <- input$steepness
    SSB2 <- c(SSB_init2,rep(0,term))
    Catch <- rep(0,term)
    N_at_age_mat <- F_tmp0 <- matrix(0,nrow = term+1,ncol=age)
    N_at_age_mat[1,] <- N_at_age
    for (y in 1:term) {
      F_tmp0[y,] <- F_at_age * exp(rnorm(1,mean = 0,sd = sd2)-sd2^2/2)
      Catch[y] <- sum(
        F_tmp0[y,]/(F_tmp0[y,]+M_at_age)*(1-exp(-F_tmp0[y,]-M_at_age))*N_at_age_mat[y,]*W_at_age
      )
      for (a in 1:age){
        if(a==1) {
          N_at_age_mat[y+1,a] <- BH2(SSB2[y],B0,R0,steepness)* exp(rnorm(1,mean = 0,sd = sd)-sd^2/2)
        }else{
          N_at_age_mat[y+1,a] <- N_at_age_mat[y,a-1] * exp(-M_at_age[a-1]-F_tmp0[y,a-1])
        }
      }
      SSB2[y+1] <- sum(N_at_age_mat[y+1,]*Q_at_age*W_at_age)
    }
    #plot(1:term,SSB2[1:term],type = "l")
    S_at_age <- c(1,rep(0,age-1))
    for (a in 1:(age-1)){
      S_at_age[a+1] <- S_at_age[a] * exp(-M_at_age[a]-F_at_age[a])
    }
    SPR1 <- sum(S_at_age * W_at_age * Q_at_age)
    #print(SPR1) # なお、a*SPR1 - b < 0で絶滅する
    #print(SPR1/SPR0) # %SPR（漁業がない状態のSPRで割ったSPR）
    #################################
    # SPRカーブの描画
    f_mulitplier <- seq(0,3,0.01)
    f_len <- length(f_mulitplier)
    SPR_res <- rep(0,f_len)
    for (i in 1:f_len){
      
      S_at_age <- c(1,rep(0,age-1))
      for (a in 1:(age-1)){
        S_at_age[a+1] <- S_at_age[a] * exp(-M_at_age[a]- f_mulitplier[i] * F_at_age[a])
      }
      SPR_res[i] <- sum(S_at_age * W_at_age * Q_at_age)/SPR0
    }
    plot(f_mulitplier,SPR_res,ylab = "%SPR",main = "%SPR",type = "l",ylim = c(0,1))
    points(x = 1,y=SPR1/SPR0,cex=2)
    
    if(input$FSPR20||input$FSPR30||input$FSPR40){#SPRシリーズを選択した場合のみ
      if (input$FSPR20){
        target_SPR <- 0.2
        F_RP <- "FSPR20"
      } else if (input$FSPR30){
        target_SPR <- 0.3
        F_RP <- "FSPR30"
      } else {
        target_SPR <- 0.4
        F_RP <- "FSPR40"
      }
      min_tmp <- min(abs(SPR_res-target_SPR))
      limit_f <- f_mulitplier[which(abs(SPR_res-target_SPR)==min_tmp)]
      points(x = limit_f,y=target_SPR,col="red",cex=2)
      abline(v=limit_f,lty=3,col="red")
    }
    #################################
    # YPR・MSYカーブの描画
    YPR <- matrix(0,nrow = f_len,ncol = age)
    MSY <- SSB_eq <- R_eq <- rep(0,f_len)
    for(i in 2:f_len){
      F_tmp_MSY <- f_mulitplier[i]*F_at_age
      S_at_age <- c(1,rep(0,age-1))
      catch_ratio <- rep(0,age)
      for (a in 1:(age-1)){
        S_at_age[a+1] <- S_at_age[a]*exp(-M_at_age[a]- F_tmp_MSY[a])
      }
      for (a in 1:age){
        catch_ratio[a] <- F_tmp_MSY[a]/(F_tmp_MSY[a]+M_at_age[a])*(1-exp(-F_tmp_MSY[a]-M_at_age[a]))
        YPR[i,a] <- catch_ratio[a]*W_at_age[a]*S_at_age[a]
      }
      SPR_tmp_MSY <- sum(S_at_age * W_at_age * Q_at_age)
      f <- function (x) BH2(x,B0,R0,steepness) - x/SPR_tmp_MSY# obtein intersection between SR−curve & RPS−line
      if (SPR_tmp_MSY * ( 4*steepness*R0/(5*steepness-1) ) - ( B0*(1-steepness)/(5*steepness-1) ) < 0 ) {#break # extinction
        SSB_eq[i] <- 0
        R_eq[i] <- 0
        MSY[i] <- 0
      } else {
        sol <- uniroot(f, c(0.00001, 1000*SSB_init2))
        SSB_eq[i] <- sol$root
        R_eq[i] <- SSB_eq[i]/SPR_tmp_MSY
        MSY[i] <- R_eq[i] * sum(YPR[i,])
      }
      
    }
    # YPR
    YPR_sum <- apply(YPR,1,sum)
    plot(f_mulitplier,YPR_sum,type = "l",ylab = "YPR",main = "YPR")
    points(x = 1, y=YPR_sum[which(f_mulitplier==1)], cex=2)
    if(input$Fmax||input$F01){#YPRシリーズを選択した場合のみ
      if (input$Fmax){
        max_tmp <- max(YPR_sum)
        limit_f <- f_mulitplier[which(YPR_sum==max_tmp)]
        points(x = limit_f,y=max_tmp,col="red",cex=2)
        abline(v=limit_f,lty=3,col="red")
        F_RP <- "Fmax"
      } else {
        YPR_slope <- (YPR_sum[2:f_len] - YPR_sum[1:(f_len-1)])
        YPR_slope_10_origin <- 0.1*YPR_slope[1]
        #cat(YPR_slope_10_origin,"\n")
        min_tmp <- min(abs(YPR_slope - YPR_slope_10_origin))
        #cat(min_tmp,"\n")
        limit_f <- f_mulitplier[which(abs(YPR_slope - YPR_slope_10_origin)==min_tmp)]
        #cat(limit_f,"\n")
        points(x = limit_f, y=YPR_sum[which(f_mulitplier==limit_f)],col="red",cex=2)
        abline(v=limit_f,lty=3,col="red")
        F_RP <- "F01"
      }
    }
    # MSY
    plot(f_mulitplier,MSY,type = "l",ylab = "sustainable yield",main = "MSY")
    points(x = 1, y=MSY[which(f_mulitplier==1)],cex=2)
    B_MSY <- SSB_eq[which(MSY==max(MSY))]
    
    if(input$Fmsy){
      limit_f <- f_mulitplier[which(MSY==max(MSY))]
      points(x = limit_f, y=max(MSY),col="red",cex=2)
      abline(v=limit_f,lty=3,col="red")
      F_RP <- "Fmsy"
    }
    #################################
    # 管理基準となる漁獲係数
    plot(age_vector,F_at_age,type = "b",xlab = "age",ylab = "F",main = "selectivity",ylim = c(0,1))
    if(!is.null(F_RP)){
      F_limit_at_age <- limit_f*F_at_age
      points(age_vector,F_limit_at_age,type = "b",col="red")
      mtext(text=F_RP,line=-2,col=2)
      management <- TRUE
    } else {
      management <- FALSE
    }
    #################################
    
    B_20 <- 0.2*B0
    B_R50 <- BH_b * R0/2 / (BH_a - R0/2) 
    
    if(input$SSB20){
      B_limit <- B_20
      B_RP <- "SSB20"
    } else if (input$R50) {
      B_limit <- B_R50
      B_RP <- "R50"
    } else if (input$Bmsy){
      B_limit <- B_MSY
      B_RP <- "Bmsy"
    } else {
      management <- FALSE
    }
    
    #################################
    # 管理
    if(management){
      SSB_init3 <- SSB2[y+1]
      term_management <- 30 
      N_at_age_mat_new <- F_tmp <- matrix(0,nrow = term_management+1,ncol=age)
      N_at_age_mat_new[1,] <- N_at_age_mat[term+1,]
      SSB3 <- c(SSB_init3,rep(0,term_management))
      Catch2 <- rep(0,term_management)
      for(y in 1:term_management){
        if(SSB3[y] < B_limit) {
          F_tmp[y,] <- F_limit_at_age * SSB3[y] / B_limit * exp(rnorm(1,mean = 0,sd = sd2)-sd2^2/2)
        }else {
          F_tmp[y,] <- F_limit_at_age * exp(rnorm(1,mean = 0,sd = sd2)-sd2^2/2)
        }
        
        Catch2[y] <- sum(
          F_tmp[y,]/(F_tmp[y]+M_at_age)*(1-exp(-F_tmp[y]-M_at_age))*N_at_age_mat_new[y,]*W_at_age
        )
        
        for (a in 1:age){
          if(a==1) {
            N_at_age_mat_new[y+1,a] <- BH2(SSB3[y],B0,R0,steepness)* exp(rnorm(1,mean = 0,sd = sd)-sd^2/2)
          }else{
            N_at_age_mat_new[y+1,a] <- N_at_age_mat_new[y,a-1] * exp(-M_at_age[a-1]-F_tmp[y,a-1])
          }
        }
        SSB3[y+1] <- sum(N_at_age_mat_new[y+1,]*Q_at_age*W_at_age)
      }
      plot(1:(term + term_management),c(SSB2[1:term],SSB3[1:term_management]),type = "l",xlab = "year",ylab = "SSB",main = "SSB_HCR",ylim = c(0,3))
      abline(v=term+1,lty=3)
      abline(h=B_limit,col="red")
      mtext(text=B_RP,line=-2,col=2)
      plot(1:(term + term_management),c(Catch[1:term],Catch2[1:term_management]),type = "l",xlab = "year",ylab = "catch",main = "catch_HCR",ylim = c(0,1.5))
      abline(v=term+1,lty=3)
      
      # kobeプロットの準備
      plot(NULL,xlim = c(0,3),ylim = c(0,3),xlab = "SSB/B_limit",ylab = "F/F_limit",main = "kobe plot")
      polygon(c(-1,1,1,-1),c(-1,-1,1,1),col="khaki1",border=NA)
      polygon(c(1,6,6,1),c(-1,-1,1,1),col="olivedrab2",border=NA)
      polygon(c(1,6,6,1),c(1,1,6,6),col="khaki1",border=NA)
      polygon(c(-1,1,1,-1),c(1,1,6,6),col="indianred1",border=NA)
      
      # 管理前をプロット
      kobe_F1 <- F_tmp0[1:term,age]/F_limit_at_age[age]
      kobe_SSB1 <- SSB2[1:term]/B_limit
      points(kobe_SSB1,kobe_F1,type = "l",col="blue")
      points(kobe_SSB1,kobe_F1,type = "p",col="blue",pch=".",cex=3)
      
      # 管理後をプロット
      kobe_F2 <- F_tmp[1:term_management,age]/F_limit_at_age[age]
      kobe_SSB2 <- SSB3[1:term_management]/B_limit
      points(kobe_SSB2,kobe_F2,type = "l")
      points(kobe_SSB2,kobe_F2,type = "p",pch=".",cex=3)
    }
    
    # 加入-産卵親魚量の関係
    curve(BH2(x,B0,R0,steepness) * exp(rnorm(1,mean = 0,sd = sd)-sd^2/2),from = 0,to = 5,ylim=c(0,1),
          xlab = "SSB",ylab = "recruitment",main="S-R relationship")
    points(SSB2[1:(term-1)],N_at_age_mat[2:term,1])
  })
  
  
}

shinyApp(ui = ui, server = server)

