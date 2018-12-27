## Kenshu 2018.12.
## 岡村
## 必要なRパッケージなど：MuMIn，TMB，Rtools
##

## データの読み込みと下調べ

dat1 <- read.csv("dat1.csv")   #  データの読み込み

head(dat1)   #  最初の6行を表示

table(dat1$count)    #  countデータの内訳を見る

plot(count~year, data=dat1)    #  年トレンドのプロット

tapply(dat1$count,dat1$year,mean)     #  countデータの年平均

tapply(dat1$count,list(dat1$plant,dat1$year), mean)         #  countデータの年/植生別の平均

plot(count~year, data=subset(dat1,plant=="tree"))      #  plant="tree"に対する年トレンド
plot(count~year, data=subset(dat1,plant=="shrub"))       #  plant="shrub"に対する年トレンド

## 正規線形回帰モデル

b <- cov(dat1$count, dat1$year)/var(dat1$year)     #  　線形回帰の傾きの最小二乗推定値
a <- mean(dat1$count)-b*mean(dat1$year)     #  　線形回帰の切片の最小二乗推定値
c(a,b)        #  　最小二乗推定値の表示

lm(count~year, data=dat1)        #    countの年変化の線形回帰      

lm(count~year*plant, data=dat1)       #    count~year+plant+year:plantの線形回帰      

# AICのシミュレーション

Res.aic <- NULL     #   結果の入れ物を作る
Sim <- 10000      #   シミュレーション回数（10000回ぐらいしないと良い結果にならない）
set.seed(1234)      #   乱数の初期値を設定
for (n in c(10,100,500,1000)){       #   サンプルサイズを10, 100, 500, 1000と変化させる
  z <- matrix(rnorm(n*Sim, 0, 1),nrow=Sim,ncol=n)      #   シミュレーション回数分データを発生（applyを使うため行列形式にしている）
  TL <- apply(z,1,function(z) integrate(function(x) dnorm(x,0,1)*dnorm(x,mean(z),sqrt(var(z)*(n-1)/n),log=TRUE),-Inf,Inf)$value)
  #  平均対数尤度   \int Q(x) P(x|\theta)．分散に(n-1)/nを掛けているのは，不偏分散を最尤推定量にするため（研修ではここ間違ってました．失礼しました）．dnormは正規分布の確率密度関数 dnorm(x, 平均, 標準偏差)                                       
  EL <- apply(z,1,function(z) mean(dnorm(z,mean(z),sqrt(var(z)*(n-1)/n),log=TRUE)))
  #  標本に対する対数尤度の平均．
  Res.aic <- cbind(Res.aic, n*(EL-TL))
  #   2つの差はパラメータ数（この場合，平均と分散の2個）になることが期待される
}
colnames(Res.aic) <- c(10,100,500,1000)     #  列の名前をサンプルサイズにする
colMeans(Res.aic)       #  Res.aicはSim × 4の行列になっていて，列平均をとっている．これの値が2になっていれば，AICの-2*(対数尤度-パラメータ数)のパラメータ数を引くというところが正しいというのが確かめられたことになる

res.n0 <- lm(count~year,data=dat1)     #  count~yearの線形回帰の結果をres.n0の中に格納
res.n1 <- lm(count~year*plant,data=dat1)        #  count~year*plantの線形回帰の結果をres.n1という変数の中に格納
logLik(res.n0);logLik(res.n1)        #  res.n0とres.n1の対数尤度を表示
AIC(res.n0,res.n1)        #  res.n0とres.n1のAICを表示

plot(count~year,data=dat1,pch=as.numeric(plant)+1,col=as.numeric(plant)+1)　　　　       # 散布図を描く．データ点の形や色をplantで変える
year <- 0:9     #  年 0~9を変数として定義しておく
lines(year,res.n1$coef[1]+res.n1$coef[2]*year,col="red")     #  linesでモデルから得られた線を散布図に上書きする（plant="shrub"のとき）
lines(year,res.n1$coef[1]+res.n1$coef[3]+(res.n1$coef[2]+res.n1$coef[4])*year,col="green",lty=2)     #  linesでモデルから得られた線を散布図に上書きする（plant="tree"のとき）

library(MuMIn)     # モデル選択を行うためのパッケージ．installしてない人は事前にinstall.packages("MuMIn")としてください．
options(na.action = "na.fail") # 欠測値の取り扱いに関する設定（これをしないとdredgeが動かない）．欠測値があったときに，欠測値を無視しない．
dredge(res.n1,rank="AIC")   #  フルモデルから変数を除いたすべての組み合わせでモデルを比較．rank="AIC"でAICを使う（defaultはAICc）

## ポアソン回帰

# ポアソン回帰の尤度関数
Poisson.reg <- function(p,dat){
  lambda <- exp(p[1]+p[2]*dat$year)      # ポアソンの平均ラムダはいつも正なので，線形モデルを指数変換する．log(lambda)=a+b*yearとなるので，log link関数という．
  -sum(dat$count*log(lambda)-lambda)     # optimなど最適化関数は最小値を求めるので，対数尤度にマイナス符号をつける
}

optim(c(log(6),-0.07),Poisson.reg,dat=dat1,method="BFGS")    #  optim(パラメータの初期値, 最適化する関数, データ, method=最適化の方法) 

glm(count~year,family=poisson,data=dat1)    #  count~yearのポアソン回帰   family=poissonでポアソン分布を指定

glm(count~year*plant,family=poisson,data=dat1)   #  count~year+plant+year:plantのポアソン回帰

res.p1 <- glm(count~year*plant,family=poisson,data=dat1)     #  フルモデルをres.p1に格納
dredge(res.p1,rank="AIC")      # AICでモデル選択

#  正規線形回帰の場合と同様，ポアソン回帰のあてはまりの図を作成
plot(count~year,data=dat1,pch=as.numeric(plant)+1,col=as.numeric(plant)+1)
year <- 0:9
lines(year,exp(res.p1$coef[1]+res.p1$coef[2]*year),col="red")
lines(year,exp(res.p1$coef[1]+res.p1$coef[3]+(res.p1$coef[2]+res.p1$coef[4])*year),col="green",lty=2)

## ロジスティック回帰

# ロジスティック回帰の対数尤度関数
Logit.reg <- function(p,dat){
  prob <- 1/(1+exp(-(p[1]+p[2]*dat$year)))     #  確率は0~1の間の数字なので，こういう変換をしてやる．逆関数はlogit=log(p/(1-p))
  -sum(dat$count*log(prob)+(10-dat$count)*log(1-prob))     # 二項分布の式から
}

optim(c(0,0),Logit.reg,dat=dat1,method="BFGS")    # optimで最適化

glm(cbind(count,10-count)~year,family=binomial,data=dat1)     # glmだと cbind(成功回数，失敗回数) ~ yearのような書き方にする．他に，0と1だけの結果ならz ~ yearのような書き方も可．

res.l1 <- glm(cbind(count,10-count)~year*plant,family=binomial,data=dat1)     # ロジスティック回帰のフルモデル
dredge(res.l1,rank="AIC")       #  AICでモデル選択

## N-混合モデル

# N混合モデルの発見関数の部分（二項分布モデル）
detect.f <- function(p, dat, N, n=10){
  y <- dat$count    #  応答変数
  X <- cbind(1,as.numeric(dat$plant)-1)     #  plantはカテゴリーなので，0 (shrub) か1 (tree)に変換してやる．model.matrix関数を使うこともできる
  r <- c(1/(1+exp(-X%*%p)))    #  1匹あたりの発見確率をplantとリンクする．X%*%pは行列とベクトルの掛算
  
  dbinom(y,n,1-(1-r)^N)     #  二項確率モデル．個体はNいるので，少なくとも1回発見する確率は1-(1-r)^Nとなる．
}

log(detect.f(c(1,1),dat1,2*dat1$count))[1:5]    #  発見モデルがちゃんと動くかテスト

# N混合モデルの個体数変動の部分（ポアソン分布モデル）
pop.f <- function(p, dat, N){
  X <- cbind(1,dat$year)    #  説明変数の行列
  lambda <- c(exp(X%*%p))    #  平均個体数（年で変化する）
  
  dpois(N,lambda)     #  ポアソン尤度
}

# N混合モデルの尤度関数
popdet.f <- function(p, dat, max.N=100, n=10){
  Pop <- function(i) {
    pop1 <- pop.f(p[3:4],dat[i,],0:max.N)       #  個々のデータに対して，ポアソン確率を0からmax.Nまで計算する
    pop1 <- pop1/sum(pop1)    # max.N=100で打ち切ったので，確率にするため全部の和で割ってやる
    pop1
  }
  like <- sapply(1:nrow(dat), function(i) sum(detect.f(p[1:2],dat[i,],0:max.N,n)*Pop(i)))      #  0~max.Nを足し合わせた周辺尤度関数を各データ（行）に対して計算
  -sum(log(like))    # 対数尤度にして足し合わせる
}

(res.nm <- optim(c(0,0,0,0),popdet.f,dat=dat1,max.N=100,method="BFGS"))     # optimで最適化  ()でくくると結果が表示される

table(dat1$plant,dat1$year)      #  plantのyearによる変化を見てやって，個体数変化の解釈をする

# N混合モデルの曲線をプロットしてやる
plot(count~year,data=dat1)      
lines(year,exp(res.nm$par[3]+res.nm$par[4]*year))

## 精度計算と信頼区間（対数正規モデル） 

# 精度計算のため，hessian=TRUEとしてヘッセ行列を取り出す
(res.nm <- optim(c(0,0,0,0),popdet.f,dat=dat1,max.N=100,method="BFGS",hessian=TRUE))

# ヘッセ行列の逆行列はパラメータ推定値の分散共分散行列になっている
(var.nm <- solve(res.nm$hessian))

# パラメータは個体数をlog scaleで見たものなので元のスケールに戻してやる
pred.N <- function(p,year) exp(p[3]+p[4]*year)
(E.N <- pred.N(res.nm$par,year))

# 個体数推定値をパラメータで微分したもの（数値微分）を計算して，デルタ法で個体数の分散や標準誤差，CVを計算                 
DN <- function(p,year,d=0.0001){
  dp <- diag(d,length(p))
  apply(dp,1,function(x) (pred.N(p+x,year)-pred.N(p-x,year))/(2*d))
}
dN <- DN(res.nm$par,0:9)
(var.N <- apply(dN, 1, function(x) x%*%var.nm%*%x))
(se.N <- sqrt(var.N))

(CV.N <- se.N/E.N)

# 上のCVを使って対数正規の95%信頼区間をプロットしてやる                
plot(count~year,data=dat1)
lines(year,E.N)
alpha <- 0.05
C.ci <- exp(qnorm(1-alpha/2)*sqrt(log(1+CV.N^2)))
Lo.lim <- E.N/C.ci
Up.lim <- E.N*C.ci
arrows(year,Lo.lim,year,Up.lim,angle=90,code=3,length=0.03,col="blue") 

## ブートストラップ信頼区間
# 上では対数正規を使ったが，個体数に対数正規の仮定をおかないで信頼区間を作りたい

# ノンパラメトリックブートストラップ
B <- 1000
b.boot <- NULL
for (i in 1:B){
  res.n0.b <- lm(count~year, data=dat1[sample(100,replace=TRUE),])    # リサンプリングしたデータに線形回帰をフィット
  b.boot <- c(b.boot, res.n0.b$coef[2])     #  傾きの推定値をとっておく
}
quantile(b.boot,probs=c(0.025,0.975))     #  quantile関数で95%信頼区間を出す

# 残差ブートストラップ
# 上のノンパラブートストラップでは，説明変数の出現頻度も変わってしまう，それが望ましくない場合もある．そこで残差をブートストラップする．
resid.n0 <- residuals(res.n0)     #  残差を計算
b.rb <- NULL
for (i in 1:B){
  count.rb <- predict(res.n0)+resid.n0[sample(100,replace=TRUE)]      #  リサンプリングした残差を予測値に足して，新たなデータを作る（ここも研修のとき間違ったコードでした．すいません）
  res.n0.rb <- lm(count.rb~year, data=dat1)     # 新たなデータに線形回帰モデルをフィット
  b.rb <- c(b.rb, res.n0.rb$coef[2])      #  傾きをとっておく
}
quantile(b.rb,probs=c(0.025,0.975))     #  95%信頼区間

#  パラメトリックブートストラップ
#  ポアソンの場合，カウントデータが応答変数になるので残差ブートストラップは望ましくない
res.p0 <- glm(count~year,family=poisson,data=dat1)     #  ポアソン回帰の結果
pred.p0 <- predict(res.p0,type="response")     # 予測値を計算してやる．predictはdefaultで線形予測子を返すので，type="response"として元のスケールに戻してやる（正規分布の場合はidentity linkなのでこの操作必要なかった）
b.pb <- NULL
for (i in 1:B){
  count.pb <- rpois(100,pred.p0)      #  予測値を使ってポアソン乱数を発生
  res.p0.pb <- glm(count.pb~year,family=poisson,data=dat1)      #  生成したデータにポアソン回帰を実行
  b.pb <- c(b.pb, res.p0.pb$coef[2])    #  傾きを記録する
}
quantile(b.pb,probs=c(0.025,0.975))     #  95%信頼区間

# N混合モデルの確率を計算する関数
dnm <- function(x, p, year, plant, max.N=100, n=10){
  lambda <- exp(p[3]+p[4]*year)
  pop1 <- dpois(0:max.N,lambda)
  pop1 <- pop1/sum(pop1)
  pop1
  r <- 1/(1+exp(-(p[1]+p[2]*plant)))
  sum(dbinom(x,n,1-(1-r)^(0:max.N))*pop1)
}

# 各データに対して，発見数が0から10までの場合の確率を計算してやる                
prob.nm <- sapply(1:100, function(i) sapply(0:10,function(x) dnm(x,res.nm$par,as.numeric(dat1$year[i]),as.numeric(dat1$plant[i])-1)))
cprob.nm <- apply(prob.nm,2,cumsum)     #  累積確率分布を計算

dat1.cb <- dat1
d <- 10^(-10)     #  数値計算の精度の調整のために必要
b.nm <- NULL
for (i in 1:B){
  z.cb <- runif(100,0,1)      # 一様乱数を発生
  count.nm <- sapply(1:100, function(x) which(cprob.nm[,x] >= z.cb[x]-d)[1]-1)     #  上で計算した累積確率との比較で，データ（0~10）を作り出してやる
  dat1.cb$count <- count.nm      #  できたデータを応答変数にする
  res.nm.cb <- optim(c(0,0,0,0),popdet.f,dat=dat1.cb,max.N=100,method="BFGS")     # N混合モデルによるパラメータ推定
  b.nm <- rbind(b.nm, exp(res.nm.cb$par[3]+res.nm.cb$par[4]*year))      # 0~9年の個体数推定値を計算
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
