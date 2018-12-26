
##### R code for CPUE standardization in Shigen-kanri kensyu, 2018

## Author: Shota Nishijima (NRIFS)

### Example 1: Shako

dat <- read.csv("shako.csv", header=TRUE)
dat <- subset(dat, Year<2018)  #2018年を除く
nrow(dat) #N = 192 = 12 months x 16 years
head(dat)

par(mfrow=c(4,4),mar=c(2,2,2,1))
for(i in unique(dat$Year)){
  plot(CPUE ~ Month, data = subset(dat,Year==i),type="l",lwd=2,
       ylim=c(0,max(dat$CPUE)),xlab="",ylab="",main=i,cex.main=2)
}

years <- unique(dat$Year)
nyear <- length(years)

dat$Year <- as.factor(dat$Year) #カテゴリカル変数に変換
dat$Month <- as.factor(dat$Month) #カテゴリカル変数に変換

## normal distribution

norm <- glm(CPUE ~ Year + Month, family = gaussian(), data = dat)
summary(norm)
as.data.frame(coef(norm))

## obtaining year trend

standardCPUE <- c()
for (i in 1:nyear) {
  if (i==1) standardCPUE[i] <- norm$coefficients[1] #1年目は切片の係数
  else standardCPUE[i] <- 
      norm$coefficients[1] + norm$coefficients[i] # 2年目以降は切片+回帰係数
}
standardCPUE <- standardCPUE/mean(standardCPUE) #scaled by mean

## plot

# nominal CPUE (mean)
nominalCPUE <- tapply(dat$CPUE, dat$Year, mean)
nominalCPUE <- nominalCPUE/mean(nominalCPUE) #scaling by mean

CPUEsummary <- cbind(nominal = nominalCPUE, Normal = standardCPUE)
dev.off()
matplot(years,CPUEsummary,ylab="CPUE",xlab="Year",
        type="b",pch=1:ncol(CPUEsummary),lwd=3,cex=2,cex.axis=1.5,
        ylim=c(0,2),cex.lab=2)

# no intercept
norm <- glm(CPUE ~ 0 + Year + Month, family = gaussian(), data =dat)
summary(norm)

standardCPUE <- norm$coefficients[1:nyear] #係数をそのまま取ればいい
standardCPUE <- standardCPUE/mean(standardCPUE)


# diagnostics

hist(resid(norm), freq=FALSE)
X <- seq(min(resid(norm)) * 1.3, max(resid(norm)) * 1.3, length = 200)
points(X, dnorm(X, 0, sqrt(mean(resid(norm)^2))), col = 2, lwd = 3, type = "l")

# test for normality
shapiro.test(resid(norm)) # not significant

par(mfrow=c(2,2), mar=c(5,5,3,2))
plot(norm,cex=1.5,lwd=3,cex.lab=2)

## lognormal distribution
lognorm <- glm(log(CPUE) ~ 0 + Year + Month, family = gaussian(), data=dat)
summary(lognorm)

# Year trend
standardCPUE <- exp(lognorm$coefficients[1:nyear]) #exp
standardCPUE <- standardCPUE/mean(standardCPUE)

# plot
CPUEsummary <- cbind(CPUEsummary, Lognormal = standardCPUE)
par(mfrow=c(1,1), mar=c(5,5,3,2))
matplot(years,CPUEsummary,ylab="CPUE",xlab="Year",
        type="b",pch=1:ncol(CPUEsummary),lwd=2,cex=2,cex.axis=1)

# check for normality
hist(resid(lognorm), freq=FALSE)
X <- seq(min(resid(lognorm)) * 1.3, max(resid(lognorm)) * 1.3, length = 200)
points(X, dnorm(X, 0, sqrt(mean(resid(lognorm)^2))), col = 2, lwd = 3, type = "l")

shapiro.test(resid(lognorm)) # significant

qqnorm(resid(lognorm), cex = 2)
qqline(resid(lognorm), lwd = 3,col=2)

# gamma distribution

gamma <- glm(CPUE ~ 0 + Year + Month, family = Gamma("log"), data = dat)

# Year trend
standardCPUE <- exp(gamma$coefficients[1:nyear]) #exp
standardCPUE <- standardCPUE/mean(standardCPUE)

# Homework; model diagnostics for the gamma model

# Comparison by AIC
AIC(gamma) # gamma
as.numeric(-2*(logLik(lognorm)-sum(log(dat$CPUE))) + 
             2*(norm$rank+1)) #AIC for lognormal

# plot
CPUEsummary <- cbind(CPUEsummary, Gamma = standardCPUE)
matplot(years,CPUEsummary,ylab="CPUE",xlab="Year",
        type="b",pch=1:ncol(CPUEsummary),lwd=2,cex=2,cex.axis=1)


# Poisson distribution

pois <- glm(Catch ~ Year + Month + offset(log(Effort)), family = poisson("log"),data = dat)
summary(pois)

# overdispersion test 
# install.packages("AER")
library(AER) #library AERの読み込み
dispersiontest(pois) #significantly overdispersed

# negative binomial
# install.packages("MASS")
library(MASS)
nb <- glm.nb(Catch ~ Year + Month + offset(log(Effort)), data = dat)
summary(nb)
AIC(pois,nb)

dat$pois.pred <- predict(pois,type="response")/dat$Effort
dat$nb.pred <- predict(nb,type="response")/dat$Effort
par(mfrow=c(4,4),mar=c(3,3,3,2))
for(i in years){
  subdat <- subset(dat,Year==i)
  matplot(cbind(subdat$CPUE,subdat$pois.pred,subdat$nb.pred),
          type=c("p","l","l"),lwd=2,pch=16,lty=1,cex.main=2,cex.axis=2,
       ylim=c(0,max(dat$CPUE)),xlab="",ylab="",main=i)
}

### interaction between Year and Month
interact0 <- glm(log(CPUE) ~ Year*Month, data = dat)
summary(interact0)
# overfitting

# Regard Month as a continuous variable
dat$Month <- as.numeric(dat$Month)
interact <- glm(log(CPUE) ~ Year * Month, data = dat)
summary(interact)

dat$lognorm.pred <- exp(predict(lognorm))
dat$interact.pred <- exp(predict(interact))

interact2 <- glm(log(CPUE) ~ Year * Month + I(Month^2), data = dat)
interact3 <- glm(log(CPUE) ~ Year * Month + Year * I(Month^2), data = dat)
AIC(lognorm, interact, interact2, interact3)


dat$interact.pred3 <- exp(predict(interact3))
par(mfrow=c(4,4),mar=c(3,3,3,2))
for(i in years){
  subdat <- subset(dat,Year==i)
  matplot(cbind(subdat$CPUE,subdat$lognorm.pred,subdat$interact.pred,subdat$interact.pred3),
          type=c("p","l","l","l"),lwd=3,pch=16,lty=1,cex.main=2,
          ylim=c(0,max(dat$CPUE)),xlab="",ylab="",main=i,cex.axis=2)
}


## obtain Year trend 

new.dat <- expand.grid(Year=unique(dat$Year), Month=unique(dat$Month))
new.dat$pred.cpue <- exp(predict(interact3, newdata=new.dat))
standardCPUE <- tapply(new.dat$pred.cpue, new.dat$Year, mean)
standardCPUE <- standardCPUE / mean(standardCPUE)

par(mfrow=c(1,1),mar=c(5,5,3,3))
matplot(years, cbind(nominalCPUE,standardCPUE),type="b",pch=c(1,19),
        cex=3,lwd=3,xlab="Year",ylab="CPUE",cex.axis=2,cex.lab=2)

### calculating bootstrap confidence intervals

## data resampling
set.seed(12345)
nsim <- 1000 #ブートストラップ回数
boot.cpue <- sapply(1:nsim, function(i){
  boot.dat <- dat[sample(1:nrow(dat),nrow(dat),replace=TRUE),] #データのリサンプリング
  boot.res <- update(interact3, data=boot.dat) #GLMの更新
  pred.cpue <- exp(predict(boot.res, newdata=new.dat)) #CPUEの予測
  staCPUE <- tapply(pred.cpue, new.dat$Year, mean) #標準化CPUEの導出
  staCPUE / mean(staCPUE)
})

ci <- c(0.025,0.975)
plot(years,standardCPUE,type="n",pch=c(19),cex=2,lwd=2,ylim=c(0,5),col="red",
     xlab="",ylab="",xaxt="n",yaxt="n")
polygon(x=c(years,rev(years)),
        y=c(apply(boot.cpue,1,quantile,probs=ci[1]),rev(apply(boot.cpue,1,quantile,probs=ci[2]))),
        border="white",col="lightblue")
par(new=T)
plot(years,standardCPUE,type="b",pch=c(19),cex=2,lwd=3,ylim=c(0,5),col="blue",
     xlab="Year",main="95% Confidence interval (data resampling)",
     cex.main=2,cex.axis=2,cex.lab=2,ylab="Standardized CPUE")

### residual resampling
set.seed(12345)
boot.cpue2 <- sapply(1:nsim, function(i){
  boot.resid <- sample(interact3$resid,nrow(dat),replace=TRUE) #残差のリサンプリング
  boot.dat <- dat
  boot.dat$CPUE <- exp(interact3$fitted.values + boot.resid) #予測値に残差を足す
  boot.res <- update(interact3, data=boot.dat) #後は同じ
  pred.cpue <- exp(predict(boot.res, newdata=new.dat))
  staCPUE <- tapply(pred.cpue, new.dat$Year, mean)
  staCPUE / mean(staCPUE)
})

plot(years,standardCPUE,type="n",pch=c(19),cex=2,lwd=2,ylim=c(0,5),col="red",
     xlab="",ylab="",xaxt="n",yaxt="n")
polygon(x=c(years,rev(years)),
        y=c(apply(boot.cpue2,1,quantile,probs=ci[1]),rev(apply(boot.cpue2,1,quantile,probs=ci[2]))),
        border="white",col="lightblue")
par(new=T)
plot(years,standardCPUE,type="b",pch=c(19),cex=2,lwd=3,ylim=c(0,5),col="blue",
     xlab="Year",main="95% Confidence interval (residual resampling)",
     cex.main=2,cex.axis=2,cex.lab=2,ylab="Standardized CPUE")

### posterior distribution of parameters

# install.packages("arm")
library(arm)　#package"arm"の読み込み
par.post <- sim(interact3, n.sims = 1000) #事後分布の計算
boot.coef <- coef(par.post)　#係数の取り出し
new.dat$CPUE <- rep(1,length(new.dat))
sim.dat <- model.matrix(interact3$formula, data = new.dat) %*% t(boot.coef)
#予測値の計算
#model.matrixはカテゴリ変数をダミー変数化したり係数にかける値を算出する関数
sim.dat <- exp(sim.dat)
boot.cpue3 <- apply(sim.dat,2, function(x) {
  staCPUE <- tapply(x, new.dat$Year, mean)
  staCPUE / mean(staCPUE)
  })

par(mfrow=c(1,3),mar=c(4,4,2,2),lwd=3)
hist(coef(par.post)[,1],main="")
hist(coef(par.post)[,2],main="")
hist(coef(par.post)[,11],main="")

par(mfrow=c(1,1))
plot(years,standardCPUE,type="n",pch=c(19),cex=2,lwd=2,ylim=c(0,5),col="red",
     xlab="",ylab="",xaxt="n",yaxt="n")
polygon(x=c(years,rev(years)),
        y=c(apply(boot.cpue3,1,quantile,probs=ci[1]),rev(apply(boot.cpue3,1,quantile,probs=ci[2]))),
        border="white",col="lightblue")
par(new=T)
plot(years,standardCPUE,type="b",pch=c(19),cex=2,lwd=2,ylim=c(0,5),col="blue",
     xlab="Year",main="95% Confidence interval (posterior distribution)",
     cex.main=2,cex.axis=2,cex.lab=2,ylab="Standardized CPUE")


### GAM 

# install.packages("mgcv")
library(mgcv)
dat$Month <- as.numeric(dat$Month)

spline <- gam(log(CPUE) ~ s(as.numeric(1:nrow(dat))), data = dat)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(spline, resid=TRUE,cex=1,pch=1,xlab="Year-month",ylab="log(CPUE)",
     cex.lab=2,cex.axis=2)

sp.ls <- exp(seq(log(0.00001), log(1), length=100))
spline.ls <- lapply(sp.ls, function(i) {
  gam(log(CPUE) ~ s(as.numeric(1:nrow(dat))), data = dat,
      sp = i)
  })

par(mfrow=c(2,2),mar=c(5,5,2,2))
plot(sp.ls, sapply(1:100, function(i) spline.ls[[i]]$gcv.ubre), log="x",
     xlab = "sp", ylab="GCV",cex.lab=2)
points(spline$sp,spline$gcv.ubre,cex=2,col=2,pch=19)
plot(spline.ls[[100]],resid=TRUE,cex=1,pch=1,xlab="Year-month",ylab="log(CPUE)",
     cex.lab=2,cex.axis=2)
plot(spline.ls[[1]],resid=TRUE,cex=1,pch=1,xlab="Year-month",ylab="log(CPUE)",
     cex.lab=2,cex.axis=2)
plot(spline,resid=TRUE,cex=1,pch=1,xlab="Year-month",ylab="log(CPUE)",
     cex.lab=2,cex.axis=2)

gam <- gam(log(CPUE) ~ Year+s(Month), data = dat) #no interaction
gam2 <- gam(log(CPUE) ~ Year + s(Month, by=Year), data = dat) # interaction 
AIC(gam,gam2)
c(gam$gcv.ubre,gam2$gcv.ubre)

dat$gam.pred <- exp(predict(gam))
dat$gam.pred2 <- exp(predict(gam2))
par(mfrow=c(4,4),mar=c(2,2,2,1),lwd=2)
for(i in years){
  subdat <- subset(dat,Year==i)
  matplot(cbind(subdat$CPUE,subdat$gam.pred,subdat$gam.pred2),
          type=c("p","l","l"),lwd=2,pch=16,lty=1,
          ylim=c(0,max(dat$CPUE)),xlab="",ylab="",main=i)
}

### GLM vs. GAM (cross validation)

set.seed(1234)
nfolds <- 10 #10-fold cross validation
dat$id <- sample(1:nfolds, nrow(dat), replace=TRUE)
rmse <- sapply(1:nfolds, function(i) {
  train.dat <- subset(dat, id != i)
  test.dat <- subset(dat, id == i)
  glm.cv <- predict(update(interact3, data=train.dat),newdata = test.dat)
  glm.rmse <- sqrt(mean((log(test.dat$CPUE) - glm.cv)^2))
  gam.cv <- predict(update(gam2, data=train.dat),newdata = test.dat)
  gam.rmse <- sqrt(mean((log(test.dat$CPUE) - gam.cv)^2))
  c(glm.rmse,gam.rmse)
  })

rmse
rowMeans(rmse)

## Homework: Calculate year trend and its confidence interval using the residual bootstrap


### Example 2: Minamimaguro

dat <- read.csv("minamimaguro.csv", header=T)
dat$TARGET_SPECIES <- NULL
dat <- na.omit(dat) #NAを除く

dat$LONGITUDE <- ifelse(dat$LONGITUDE < -70, dat$LONGITUDE+360, dat$LONGITUDE)

dat$Catch <- as.integer(dat$NUMBER_OF_SBT_RETAINED)
dat$Effort <- dat$NUMBER_OF_HOOKS
dat$CPUE <- dat$NUMBER_OF_SBT_RETAINED / dat$NUMBER_OF_HOOKS

dat <- subset(dat, YEAR >= 1969 & 
                CCSBT_STATISTICAL_AREA %in% as.factor(4:9) & 
                MONTH >= 4 & MONTH <= 9)
#データを限定

dat$Year <- as.factor(dat$YEAR)
dat$StatArea <- dat$CCSBT_STATISTICAL_AREA
dat$Month <- as.factor(dat$MONTH)
dat$Area <- dat$StatArea
dat$Country <- as.factor(dat$COUNTRY_CODE)

dat.p <- subset(dat, CPUE>0)

1-nrow(dat.p)/nrow(dat)

### zero-inflated negative binomial model

library(pscl) #パッケージの読み込み
mod.zinb <- zeroinfl(Catch ~ Year+Month+Area+offset(log(Effort))|Year+Month+Area,
                     data = dat, dist="negbin")
cbind(mod.zinb$coefficients$count,mod.zinb$coefficients$zero)

par(mfrow=c(1,1),lwd=3,mar=c(5,5,5,2))
library(MASS)
library(DHARMa)
mod.nb <- glm.nb(Catch ~ Year + Month + Area + offset(log(Effort)), 
                 data = dat) # 負の二項分布
simresid.nb <- simulateResiduals(mod.nb) #残差をsimulate
testZeroInflation(simresid.nb) #0となる

### delta GLM

binom.model <- glm(CPUE>0 ~ Year + Area + Month, data = dat,
                   family=binomial("logit"))
gamma.model <- glm(CPUE ~ Year + Area + Month, data = dat.p, 
                   family=Gamma("log"))

# year trend
nominal <- tapply(dat$CPUE, dat$Year, mean)
nominal <- nominal / mean(nominal)

new.dat <- expand.grid(Year=unique(dat$Year), 
                       Area=unique(dat$Area), 
                       Month=unique(dat$Month))

new.dat$p <- predict(binom.model, newdata=new.dat, type = "response")
new.dat$mu <- predict(gamma.model, newdata=new.dat, type = "response")
new.dat$y <- new.dat$p * new.dat$mu 
standard <- tapply(new.dat$y, new.dat$Year, mean)
standard <- standard / mean(standard)
#本来は面積重みづけしたほうが良いと思います

par(mar=c(5,5,3,2))
matplot(unique(dat$Year),cbind(nominal,standard), type="l",lwd=3,
        ylab="CPUE",xlab="",cex.lab=2,cex.axis=2)

# Diagnostics for binomial model
simresid.b <- simulateResiduals(binom.model)
plot(simresid.b)


### delta-GLM-tree

source("tree-deltaGLM1.9.r")

res.tree <- tree.delta.glm(binom.model, gamma.model, criteria = "AIC", trend = FALSE, lat.name = "LATITUDE", lon.name = "LONGITUDE",
                           factor = c("Year", "Area", "Month"), min.lon = -20, max.lon = 190, min.lat = -50, max.lat = -30,
                           seq.lon = 10, seq.lat = 5, max.narea = 100, parallel = FALSE)

trend <- update.delta.glm(res.tree$best$binom, res.tree$best$posi, trend = TRUE,
                          factor = c("Year", "Area", "Month"))


#year trend
par(mfrow=c(1,1),lwd=3)
matplot(unique(dat$Year),cbind(nominal,trend$y.trend), type="l",lwd=3,
        ylab="CPUE",xlab="",cex.lab=2,cex.axis=2)


# mapping
library(maps)

dat$LONGITUDE2 <- ifelse(dat$LONGITUDE<0, dat$LONGITUDE+360,dat$LONGITUDE)

par(lwd=1)
plot(NA,xlim=c(0,360),ylim=c(-60,60),
     cex=0.8,col=1,xaxs="i",yaxs="i",xaxt="n",yaxt="n",ann=F)
map('world2',add=T,fill=T,col=gray(0.8))
par(new=T)
plot(cbind(dat$LONGITUDE2,dat$LATITUDE),
     cex=2,pch=19,col=res.tree$area[[2]][,3]+1,
     xlab="Longitude",ylab="Latitude",
     xlim=c(0,360),ylim=c(-60,60))


### Example 3: Sawara

dat0 <- read.csv("data.sawara.csv", header=TRUE)
head(dat0)
nrow(dat0) # N = 48613

dat0 <- dat0[order(dat0$year,dat0$month,dat0$day),] #並び替え

dat <- subset(dat0, gyokind == "定置大型") #大型のみ
dat <- subset(dat, meigara != "サワラ銘柄不明") #
nrow(dat) #N = 30128

dat$meigara <- ifelse(dat$meigara=="サワラ小","sagoshi","sawara")
dat$year <- as.factor(dat$year)
dat$month <- as.factor(dat$month)
dat$area <- sapply(1:nrow(dat), function (i) which(sort(unique(dat$netnumber))==dat$netnumber[i]))

dat.s <- subset(dat, meigara=="sagoshi") # for small size
dat.l <- subset(dat, meigara=="sawara") # for large size
nrow(dat.s) #N = 17573
nrow(dat.l) #N = 12555

dat.l$gyokind <- dat.l$gyoshu <- dat.l$meigara <- dat.l$seqno <- NULL
head(dat.l)

## GLMM

# install.packages("lme4")
library(lme4) #lme4パッケージの読み込み
glmm.res <- lmer(log(catch) ~ year + month + (1|area), data=dat.l)
summary(glmm.res)


## CAR

# install.packages("CARBayes")
library(CARBayes)

n.area <- length(unique(dat.l$area))
W <- matrix(0, ncol = n.area, nrow = n.area)
for (i in 1:(n.area-1)) W[i,i+1] <- W[i+1,i] <- 1
car.bayes <- S.CARmultilevel(formula = log(catch) ~ year + month, data = dat.l,
                             family = "gaussian", #glmのように式を指定
                             W = W, ind.area = as.numeric(dat.l$area),
                             #近傍間の重みの関係を表す行列と各データの地点を指定
                             burnin=2000, n.sample=5000, thin=3)
car.bayes$summary.results

if (0) { #複数のMCMC samplingによる収束判定）
  car.bayes.list <- lapply(1:3, function(i) {
    S.CARmultilevel(formula = log(catch) ~ year + month, data = dat.l,
                    family = "gaussian",
                    W = W, ind.area = as.numeric(dat.l$area),
                    burnin=2000, n.sample=5000, thin=3)
  })
  
  Rhat <- function(samples){
    term2 <-var(apply(samples, 2, mean))
    varW <-mean(apply(samples, 2, var))
    n.r <-nrow(samples)
    var1 <-((n.r-1)/n.r)*varW+term2
    sqrt(var1/varW)
  }
  
  # Rhat for beta
  sapply(1:dim(car.bayes$samples$beta)[2], function(i) {
    beta.samples <- sapply(car.bayes.list, function(x) x$samples$beta[,i])
    Rhat(beta.samples)
  })
  
  # Rhat for rho
  rho.samples <- sapply(car.bayes.list, function(x) x$samples$rho)
  Rhat(rho.samples)
  
  i <- 1 #intercept
  beta.samples <- sapply(car.bayes.list, function(x) x$samples$beta[,i])
  par(mfrow=c(2,1))
  matplot(beta.samples, type = "l", main="Intercept")
  legend("topright",legend=c(paste("Rhat=",sprintf("%5.2f",Rhat(beta.samples)))),
         cex=1.2)
  matplot(rho.samples, type = "l",main="rho")
  legend("topright",legend=c(paste("Rhat=",sprintf("%5.2f",Rhat(rho.samples)))),
         cex=1.2)
}


# fixed at rho=0
car.bayes.rho0 <- S.CARmultilevel(formula = log(catch) ~ year + month, 
                                  family = "gaussian", data = dat.l, W = W, 
                                  ind.area = as.numeric(dat.l$area),rho = 0,
                                  burnin=2000, n.sample=5000, thin=3)
# fixed at rho=1
car.bayes.rho1 <- S.CARmultilevel(formula = log(catch) ~ year + month, data = dat.l, 
                                  W = W, ind.area = as.numeric(dat.l$area),rho = 1,
                                  family = "gaussian", burnin=2000, n.sample=5000, thin=3)


c(rho.est=car.bayes$modelfit["WAIC"],rho.0=car.bayes.rho0$modelfit["WAIC"],rho.1=car.bayes.rho1$modelfit["WAIC"])

par(mfrow=c(1,1))
hist(car.bayes$samples$rho)

### estiamte autocorreltated random effects by years

area.year <- expand.grid(area=sort(unique(dat.l$area)), year=sort(unique(dat.l$year)))
area.year$id <- 1:nrow(area.year)

dat.l$area.year <- sapply(1:nrow(dat.l), function(i) {
  area.year$id[area.year$year==dat.l$year[i]&area.year$area==dat.l$area[i]]
})

n.year <- length(unique(dat.l$year))
W2 <- matrix(0, ncol = n.area*n.year, nrow = n.area*n.year)
for(i in 1:n.year) {
  W2[((i-1)*n.area+1):(i*n.area),((i-1)*n.area+1):(i*n.area)] <- W
}

car.bayes2 <- S.CARmultilevel(formula = log(catch) ~ year + month, data = dat.l, 
                              W = W2, ind.area = as.numeric(dat.l$area.year),
                              family = "gaussian", burnin=2000, n.sample=5000, thin=3)
car.bayes2$summary.results
hist(car.bayes2$samples$rho)
abline(v=car.bayes2$summary.results["rho",1],col=2)

car.bayes2.rho0 <- S.CARmultilevel(formula = log(catch) ~ year + month, data = dat.l, 
                                   W = W2, ind.area = as.numeric(dat.l$area.year), rho=0,
                                   family = "gaussian", burnin=2000, n.sample=5000, thin=3)

t(data.frame(rho.est=car.bayes$modelfit["WAIC"],rho.0=car.bayes.rho0$modelfit["WAIC"],
             rho.est2=car.bayes2$modelfit["WAIC"],rho.est2.rho0=car.bayes2.rho0$modelfit["WAIC"]))

new.dat <- expand.grid(year=unique(dat.l$year),
                       month=unique(dat.l$month),
                       area=unique(dat.l$area))

coef.mat <- model.matrix(~year+month, data=new.dat)

new.dat$area.year <- sapply(1:nrow(new.dat), function(i) {
  area.year$id[area.year$year==new.dat$year[i]&area.year$area==new.dat$area[i]]
})

pred0 <- car.bayes2$samples$beta %*% t(coef.mat)
re <- sapply(1:nrow(new.dat), function(i) car.bayes2$samples$phi[,new.dat$area.year[i]])
pred <- exp(pred0 + re)

new.dat$pred.med <- apply(pred,2,median) #予測値 (median)

years <- unique(dat.l$year)
par(mfcol=c(4,4), mar=c(3,3,1,1))
for(i in years) {
  subdat <- subset(dat.l,year==i)
  boxplot(catch ~ area, dat=subdat,main=i, log = "y", ylim=range(dat.l$catch),
          xlim=c(1,16),cex.main=2)
  sub.newdat <- subset(new.dat,year==i)
  par(new=T)
  plot(tapply(sub.newdat$pred.med, sub.newdat$area, mean),log="y",ylim=range(dat.l$catch),lwd=2,type="l",col=2,
       xlim=c(1,16))
}

boot.cpue <- sapply(1:nrow(pred), function(i){
  new.dat$pred <- pred[i,]
  staCPUE <- tapply(new.dat$pred, new.dat$year, mean)
  staCPUE/mean(staCPUE)
})

nominal <- tapply(dat.l$catch, dat.l$year, mean)
nominal <- nominal / mean(nominal)

standardized <- rowMeans(boot.cpue)
standardized <- standardized / mean(standardized)

standardized.low <- apply(boot.cpue,1,function(x)quantile(x, prob=0.025))
standardized.up <- apply(boot.cpue,1,function(x)quantile(x, prob=1-0.025))

years <- as.character(unique(dat.l$year))

par(mfrow=c(1,1),mar=c(5,5,3,2),cex.lab=2,cex.axis=2)
matplot(years,cbind(nominal, standardized), type = "n", pch=c(1,16),
        ylim=c(0,3.3),xlab="Year",ylab="Scaled CPUE")
polygon(x=c(years,rev(years)),
        y=c(standardized.low,rev(standardized.up)),
        border="white",col="lightpink")
par(new=T)
matplot(years,cbind(nominal, standardized), type = "b", pch=c(1,16),
        lwd=3,cex=2,ann=T,xlab="",ylab="",ylim=c(0,3.3))

## Howework: Apply the CAR model to the data of small-sized Sawara

### fin