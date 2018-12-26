## Kenshu 2018.12.
## 岡村
## 必要なRパッケージなど：MuMIn，TMB，Rtools
##

## データの読み込みと下調べ

dat1 <- read.csv("dat1.csv")

head(dat1)

table(dat1$count)

plot(count~year, data=dat1)

tapply(dat1$count,dat1$year,mean)

tapply(dat1$count,list(dat1$plant,dat1$year), mean)

plot(count~year, data=subset(dat1,plant=="tree"))
plot(count~year, data=subset(dat1,plant=="shrub"))

## 正規線形回帰モデル

b <- cov(dat1$count, dat1$year)/var(dat1$year)
a <- mean(dat1$count)-b*mean(dat1$year)
c(a,b)

lm(count~year, data=dat1)

lm(count~year*plant, data=dat1)

# AICのシミュレーション

Res.aic <- NULL
Sim <- 10000      
set.seed(1234)
for (n in c(10,100,500,1000)){
  z <- matrix(rnorm(n*Sim, 0, 1),nrow=Sim,ncol=n)
  TL <- apply(z,1,function(z) integrate(function(x) dnorm(x,0,1)*dnorm(x,mean(z),sqrt(var(z)*(n-1)/n),log=TRUE),-Inf,Inf)$value)
  EL <- apply(z,1,function(z) mean(dnorm(z,mean(z),sqrt(var(z)*(n-1)/n),log=TRUE)))
  Res.aic <- cbind(Res.aic, n*(EL-TL))
}
colnames(Res.aic) <- c(10,100,500,1000)
colMeans(Res.aic)

res.n0 <- lm(count~year,data=dat1)
res.n1 <- lm(count~year*plant,data=dat1)
logLik(res.n0);logLik(res.n1)
AIC(res.n0,res.n1)

plot(count~year,data=dat1,pch=as.numeric(plant)+1,col=as.numeric(plant)+1)
year <- 0:9
lines(year,res.n1$coef[1]+res.n1$coef[2]*year,col="red")
lines(year,res.n1$coef[1]+res.n1$coef[3]+(res.n1$coef[2]+res.n1$coef[4])*year,col="green",lty=2)

library(MuMIn)
options(na.action = "na.fail") # 欠測値の取り扱いに関する設定（これをしないとdredgeが動かない）
dredge(res.n1,rank="AIC")

## ポアソン回帰

Poisson.reg <- function(p,dat){
  lambda <- exp(p[1]+p[2]*dat$year)
  -sum(dat$count*log(lambda)-lambda)
}

optim(c(log(6),-0.07),Poisson.reg,dat=dat1,method="BFGS") 

glm(count~year,family=poisson,data=dat1)

glm(count~year*plant,family=poisson,data=dat1)

res.p1 <- glm(count~year*plant,family=poisson,data=dat1)
dredge(res.p1,rank="AIC")

plot(count~year,data=dat1,pch=as.numeric(plant)+1,col=as.numeric(plant)+1)
year <- 0:9
lines(year,exp(res.p1$coef[1]+res.p1$coef[2]*year),col="red")
lines(year,exp(res.p1$coef[1]+res.p1$coef[3]+(res.p1$coef[2]+res.p1$coef[4])*year),col="green",lty=2)

## ロジスティック回帰

Logit.reg <- function(p,dat){
  prob <- 1/(1+exp(-(p[1]+p[2]*dat$year)))
  -sum(dat$count*log(prob)+(10-dat$count)*log(1-prob))
}

optim(c(0,0),Logit.reg,dat=dat1,method="BFGS") 

glm(cbind(count,10-count)~year,family=binomial,data=dat1)

res.l1 <- glm(cbind(count,10-count)~year*plant,family=binomial,data=dat1)
dredge(res.l1,rank="AIC")

## N-混合モデル

detect.f <- function(p, dat, N, n=10){
  y <- dat$count
  X <- cbind(1,as.numeric(dat$plant)-1)
  r <- c(1/(1+exp(-X%*%p)))
  
  dbinom(y,n,1-(1-r)^N)
}

log(detect.f(c(1,1),dat1,2*dat1$count))[1:5]

pop.f <- function(p, dat, N){
  X <- cbind(1,dat$year)
  lambda <- c(exp(X%*%p))
  
  dpois(N,lambda)
}

popdet.f <- function(p, dat, max.N=100, n=10){
  Pop <- function(i) {
    pop1 <- pop.f(p[3:4],dat[i,],0:max.N)
    pop1 <- pop1/sum(pop1)
    pop1
  }
  like <- sapply(1:nrow(dat), function(i) sum(detect.f(p[1:2],dat[i,],0:max.N,n)*Pop(i)))
  -sum(log(like))
}

(res.nm <- optim(c(0,0,0,0),popdet.f,dat=dat1,max.N=100,method="BFGS"))

table(dat1$plant,dat1$year)

plot(count~year,data=dat1)
lines(year,exp(res.nm$par[3]+res.nm$par[4]*year))

## 精度計算と信頼区間（対数正規モデル） 

(res.nm <- optim(c(0,0,0,0),popdet.f,dat=dat1,max.N=100,method="BFGS",hessian=TRUE))

(var.nm <- solve(res.nm$hessian))

pred.N <- function(p,year) exp(p[3]+p[4]*year)
(E.N <- pred.N(res.nm$par,year))

DN <- function(p,year,d=0.0001){
  dp <- diag(d,length(p))
  apply(dp,1,function(x) (pred.N(p+x,year)-pred.N(p-x,year))/(2*d))
}
dN <- DN(res.nm$par,0:9)
(var.N <- apply(dN, 1, function(x) x%*%var.nm%*%x))
(se.N <- sqrt(var.N))

(CV.N <- se.N/E.N)

plot(count~year,data=dat1)
lines(year,E.N)
alpha <- 0.05
C.ci <- exp(qnorm(1-alpha/2)*sqrt(log(1+CV.N^2)))
Lo.lim <- E.N/C.ci
Up.lim <- E.N*C.ci
arrows(year,Lo.lim,year,Up.lim,angle=90,code=3,length=0.03,col="blue") 

## ブートストラップ信頼区間

B <- 1000
b.boot <- NULL
for (i in 1:B){
  res.n0.b <- lm(count~year, data=dat1[sample(100,replace=TRUE),])
  b.boot <- c(b.boot, res.n0.b$coef[2])
}
quantile(b.boot,probs=c(0.025,0.975))

resid.n0 <- residuals(res.n0)
b.rb <- NULL
for (i in 1:B){
  count.rb <- dat1$count+resid.n0[sample(100,replace=TRUE)]
  res.n0.rb <- lm(count.rb~year, data=dat1)
  b.rb <- c(b.rb, res.n0.rb$coef[2])
}
quantile(b.rb,probs=c(0.025,0.975))

res.p0 <- glm(count~year,family=poisson,data=dat1)
pred.p0 <- predict(res.p0,type="response")
b.pb <- NULL
for (i in 1:B){
  count.pb <- rpois(100,pred.p0)
  res.p0.pb <- glm(count.pb~year,family=poisson,data=dat1)
  b.pb <- c(b.pb, res.p0.pb$coef[2])
}
quantile(b.pb,probs=c(0.025,0.975))

dnm <- function(x, p, year, plant, max.N=100, n=10){
  lambda <- exp(p[3]+p[4]*year)
  pop1 <- dpois(0:max.N,lambda)
  pop1 <- pop1/sum(pop1)
  pop1
  r <- 1/(1+exp(-(p[1]+p[2]*plant)))
  sum(dbinom(x,n,1-(1-r)^(0:max.N))*pop1)
}

prob.nm <- sapply(1:100, function(i) sapply(0:10,function(x) dnm(x,res.nm$par,as.numeric(dat1$year[i]),as.numeric(dat1$plant[i])-1)))
cprob.nm <- apply(prob.nm,2,cumsum)

dat1.cb <- dat1
d <- 10^(-10)
b.nm <- NULL
for (i in 1:B){
  z.cb <- runif(100,0,1)
  count.nm <- sapply(1:100, function(x) which(cprob.nm[,x] >= z.cb[x]-d)[1]-1)
  dat1.cb$count <- count.nm
  res.nm.cb <- optim(c(0,0,0,0),popdet.f,dat=dat1.cb,max.N=100,method="BFGS")
  b.nm <- rbind(b.nm, exp(res.nm.cb$par[3]+res.nm.cb$par[4]*year))
}
CI.cb <- apply(b.nm,2,quantile,probs=c(0.025,0.975))

plot(count~year,data=dat1)
lines(year,E.N)
arrows(year,CI.cb[1,],year,CI.cb[2,],angle=90,code=3,length=0.03,col="blue")

## Template Model Builder実例

sink("nm.cpp")
cat("// N-mixture model

#include <TMB.hpp>

/* Parameter estimation */
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data Section //
  DATA_INTEGER(Nmax);
  DATA_SCALAR(n);
  DATA_MATRIX(DAT);

  // Parameter Section //
  PARAMETER_VECTOR(P);
  
  int N=DAT.rows();
  
  matrix<Type> Pois(Nmax+1,N);
  vector<Type> D(N);
  vector<Type> Log_Lambda(N);
  vector<Type> R(N);
  vector<Type> BP(N);
  
  Type d=1e-15;
  
  Type f=0;
    
  for (int j=0;j<N;j++){  
    R(j) = 1/(1+exp(-(P(0)+P(1)*DAT(j,2))));
    Log_Lambda(j) = P(2)+P(3)*DAT(j,1);
    for (int i=0;i<=Nmax;i++){
      Pois(i,j) = exp(Type(i)*Log_Lambda(j)-exp(Log_Lambda(j))-lgamma(Type(i)+1));
      D(j) += Pois(i,j);
    }
    
    for (int i=0;i<=Nmax;i++){
      BP(j) += exp(lgamma(n+1) - lgamma(DAT(j,0)+1) - lgamma(n-DAT(j,0)+1)+DAT(j,0)*log(1-pow(1-R(j),i)+d)+(n-DAT(j,0))*log(pow(1-R(j),i)+d))*Pois(i,j)/D(j);
    }
    f -= log(BP(j));
  }
  
  return f;
}",fill=TRUE)
sink()

library(TMB)
compile("nm.cpp")
dyn.load(dynlib("nm"))

dat2 <- dat1[,-3]
dat2[,3] <- as.numeric(dat2[,3])-1
dat <- list(Nmax=100, n=10, DAT=as.matrix(dat2))
parms <- list(P=c(0,0,0,0))
obj <- MakeADFun(data=dat,parameters=parms,DLL="nm")
(res.nm.tmb <- optim(obj$par,obj$fn,obj$gr,method="BFGS",hessian=TRUE))

dat1.cb.tmb <- dat1
d <- 10^(-10)
b.nm.tmb <- NULL
for (i in 1:B){
  z.cb.tmb <- runif(100,0,1)
  count.nm.tmb <- sapply(1:100, function(x) which(cprob.nm[,x] >= z.cb.tmb[x]-d)[1]-1)
  dat2[,1] <- count.nm.tmb
  dat <- list(Nmax=100, n=10, DAT=as.matrix(dat2))
  parms <- list(P=c(0,0,0,0))
  obj <- MakeADFun(data=dat,parameters=parms,DLL="nm")
  res.nm.cb.tmb <- optim(obj$par,obj$fn,obj$gr,method="BFGS")
  b.nm.tmb <- rbind(b.nm.tmb, exp(res.nm.cb.tmb$par[3]+res.nm.cb.tmb$par[4]*year))
}
CI.cb.tmb <- apply(b.nm.tmb,2,quantile,probs=c(0.025,0.975))

plot(count~year,data=dat1)
lines(year,E.N)
arrows(year,CI.cb.tmb[1,],year,CI.cb.tmb[2,],angle=90,code=3,length=0.03,col="blue")
