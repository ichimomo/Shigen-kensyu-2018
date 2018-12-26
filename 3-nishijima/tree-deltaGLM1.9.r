##### R code for applying GLM-tree to delta-GLM ###

### Reanalyzing delta-GLM by excluding categories whose all samples have CPUE > 0 or CPUE = 0 (because of not being able to use a binomial model) 
### Calculating Yearly trend
update.delta.glm <- function(
  binom.model,
  posi.model,
  trend = FALSE, #年トレンドを計算するか
  factor = c("Year", "Area"), #年トレンドを算出する際にexpand.gridをするfactorial変数
          #Yearを最初に、Areaを二番目にする（Area-weighted averageにする場合）
  numeric = NULL, #年トレンドを算出する際にexpand.gridをする連続変数
  seq.length = 10, #expand.gridをする連続変数を何個に分けるか
  offset = NULL, #offset項に使った変数（"Effort"など）※logはいらない
  area.w = TRUE, #各海区の重みづけ（FALSEですべて同じ、TRUEで緯度経度で面積重みづけ）
  lon.name = "Lon", #データセットの経度の列名
  lat.name = "Lat", #データセットの緯度の列名
  min.lon = 140, #最小経度（データが含まれれば適当で大丈夫）
  max.lon = 185, #最大経度（データが含まれれば適当で大丈夫）
  min.lat = 30, #最小経度（データが含まれれば適当で大丈夫）
  max.lat = 50, #最大経度（データが含まれれば適当で大丈夫）
  seq.lon = 2.5, #resolution of longitude
  seq.lat = 2.5, #resolution of latitude
  scale = "mean", #指標値をscalingするのに使う関数（他には"sd"など）
  remove = FALSE, #すべて0or1のカテゴリをbinomial modelから除いてupdateするか
  data = binom.model$data,
  pull.LL = 0 #正規分布⇒対数正規分布に変換する際に対数尤度に引く値
){
  # argname <- ls()
  # arglist <- lapply(argname,function(x) eval(parse(text=x)))
  # names(arglist) <- argname
  # 
  res <- list()
  # res$input <- arglist
  
  dat.all <- data
  dat.all2 <- data.frame(dat.all, predicted = binom.model$linear.predictors)
  
  binom.model2 <- binom.model
  posi.model2 <- posi.model
  
  if (remove) {
    dat.all2 <- subset(dat.all2, abs(predicted) < 10) #全て0 or 1のfactorを除く
    binom.model2 <- try(update(binom.model, data = dat.all2), silent=TRUE)
  }

  n <- nrow(dat.all)
  if((class(binom.model2)[1] == "try-error")|| (class(posi.model2)[1] == "try-error")){
    k.binom <- k.posi <- k.sum <- NULL
    loglik.binom <- loglik.posi <- loglik.sum <- -Inf
    aic <- aicc <- bic <- Inf
  } else {
    k.binom <- binom.model$rank
    k.posi <- posi.model2$rank + 1
    k.sum <- k.binom + k.posi
    loglik.binom <- as.numeric(logLik(binom.model2))
    loglik.posi <- as.numeric(logLik(posi.model2))
    if (posi.model2$family$family == "gaussian") {
      if (!is.null(posi.model2$linear.predictors)) pull.LL <- sum(posi.model2$linear.predictors)
    }
    loglik.posi <- loglik.posi - pull.LL
    loglik.sum <- loglik.binom + loglik.posi
    aic <- -2 * loglik.sum + 2 * k.sum
    aicc <- aic + 2 * k.sum * (k.sum + 1)/(n - k.sum - 1)
    bic <- -2 * loglik.sum + k.sum * log(n)
  }
  
  res$n <- n
  res$area.name <- sort(unique(dat.all2[,factor[2]]))
  res$n.area <- length(res$area.name)
  res$binom <- binom.model2
  res$posi <- posi.model2
  res$k.binom <- k.binom
  res$k.posi <- k.posi
  res$loglik.binom <- loglik.binom
  res$loglik.posi <- loglik.posi
  res$loglik.sum <- loglik.sum
  res$AIC <- aic
  res$AICc <- aicc
  res$BIC <- bic
  
  if(isTRUE(trend)){
    grid.list <- list()
    for(i in 1:length(factor)) grid.list[[i]] <- unique(dat.all[, factor[i]])
    if (!is.null(numeric)) for(i in 1:length(numeric)) grid.list[[length(factor) + i]] <- seq(min(dat.all[, numeric[i]]), max(dat.all[, numeric[i]]), length = seq.length)
    new.dat <- expand.grid(grid.list)
    if (is.null(numeric)) {
      names(new.dat) <- factor
    } else {
      names(new.dat) <- c(factor, numeric)
    }
    
    if (!is.null(offset)) {
      new.dat <- cbind(new.dat, rep(1, nrow(new.dat)))
      names(new.dat) <- c(factor, numeric, offset)
    }
    
    new.dat$pred.logitp <- predict(binom.model, newdata = new.dat)
    new.dat$pred.p <- ifelse(new.dat$pred.logitp < -10, 0, 
           ifelse(new.dat$pred.logitp > 10, 1, 1/(1+exp(-new.dat$pred.logitp))))
    new.dat$pred.m <- ifelse(new.dat$pred.p == 0, 0, 1)
    for(i in 1:nrow(new.dat)) {
      if(new.dat$pred.m[i] == 1) {
        pred.log <- try(predict(posi.model, newdata = new.dat[i,]), silent = TRUE)
        new.dat$pred.m[i] <- ifelse(class(pred.log) == "try-error", as.numeric(0), exp(pred.log))
      }
    }
    
    new.dat$pred.y <- new.dat$pred.p * new.dat$pred.m
    
    #各年・各エリアごとの予測値(p:probability, m:positive catch rate, y:cpue=p*m) 
    p.trend.full <- tapply(new.dat$pred.p, list(new.dat$Year, new.dat$Area), mean)
    m.trend.full <- tapply(new.dat$pred.m, list(new.dat$Year, new.dat$Area), mean)
    y.trend.full <- tapply(new.dat$pred.y, list(new.dat$Year, new.dat$Area), mean)

    ### 年トレンドの算出
    seq.lon0 <- seq(min.lon, max.lon, by = seq.lon)
    seq.lat0 <- seq(min.lat, max.lat, by = seq.lat)
    
    area.vec <- sort(unique(dat.all$Area))
    if(!isTRUE(area.w)) area.w2 <- rep(1, length(area.vec)) else {
      area.w2 <- sapply(1:length(area.vec), function(j){
        lon.dis <- min(seq.lon0[seq.lon0 > max(dat.all[dat.all$Area == area.vec[j], lon.name])])-
          max(seq.lon0[seq.lon0 < min(dat.all[dat.all$Area == area.vec[j], lon.name])])
        lat.dis <- min(seq.lat0[seq.lat0 > max(dat.all[dat.all$Area == area.vec[j], lat.name])])-
          max(seq.lat0[seq.lat0 < min(dat.all[dat.all$Area == area.vec[j], lat.name])])
        lon.dis*lat.dis
      })
    } # weighted by grid sizes

    area.w3 <- area.w2/sum(area.w2)
    
    #scalingする前
    p.trend0 <- p.trend <- rowSums(p.trend.full %*% area.w3)
    m.trend0 <- rowSums(m.trend.full %*% area.w3)
    y.trend0 <- rowSums(y.trend.full %*% area.w3)
    
    #scaling
　　m.trend <- m.trend0 / get(scale)(m.trend0,na.rm=TRUE)
    y.trend <- y.trend0 / get(scale)(y.trend0,na.rm=TRUE)
    
    res$grid.data <- new.dat
    res$p.trend.full <- p.trend.full
    res$m.trend.full <- m.trend.full
    res$y.trend.full <- y.trend.full
    res$p.trend0 <- p.trend0
    res$m.trend0 <- m.trend0
    res$y.trend0 <- y.trend0
    res$area.w <- area.w3
    res$p.trend <- p.trend
    res$m.trend <- m.trend
    res$y.trend <- y.trend
  }
  
  res
}


### Operaitng area stratification by GLM-tree

tree.delta.glm <- function(
  binom.model, 
  posi.model, 
  criteria = "AIC", #モデル選択の基準（他にBIC/AICcが選択可）
  max.narea = 100, #maximum number of areas 
  trend = FALSE, #年トレンドを計算するか
  delta = 10, #delta AIC (or BIC/AICc) がこの値よりも大きくなるまで海区分けを実行
  factor = c("Year", "Area"), #年トレンドを算出する際にexpand.gridをするfactorial変数
  #Yearを最初に、Areaを二番目にする（Area-weighted averageにする場合）
  numeric = NULL, #年トレンドを算出する際にexpand.gridをする連続変数
  seq.length = 10, #expand.gridをする連続変数を何個に分けるか
  offset = NULL, #offset項に使った変数（"Effort"など）※logはいらない
  area.w = TRUE, #各海区の重みづけ（FALSEですべて同じ、TRUEで緯度経度で面積重みづけ）
  lon.name = "Lon", #データセットの経度の列名
  lat.name = "Lat", #データセットの緯度の列名
  min.lon = 140, #最小経度（データが含まれれば適当で大丈夫）
  max.lon = 185, #最大経度（データが含まれれば適当で大丈夫）
  min.lat = 30, #最小経度（データが含まれれば適当で大丈夫）
  max.lat = 50, #最大経度（データが含まれれば適当で大丈夫）
  seq.lon = 2.5, #resolution of longitude
  seq.lat = 2.5, #resolution of latitude
  scale = "mean", #指標値をscalingするのに使う関数（"sd"でも可）
  remove = FALSE, #すべて0or1のカテゴリをbinomial modelから除いてupdateするか
  cat = TRUE,
  model.ave = TRUE, #model averagingするかどうか
  parallel = FALSE, #TRUEでforeachを使ったparallel計算
  export.package = NULL, #parallel 計算するときにexportするpackage("speedglm", "glm2"など)
  cores = ifelse(parallel,parallel::detectCores(),1) #並列計算に使用するコア数
  ){
  argname <- ls()  
  arglist <- lapply(argname,function(x) eval(parse(text=x)))
  names(arglist) <- argname
  
  res <- list()
  res$input <- arglist
  
  if (parallel) {
    require(doParallel)
    type <- if (exists("mcfork", mode = "function")) "FORK" else "PSOCK"
  }
  if (posi.model$family$family == "gaussian") {
    pull.LL <- as.numeric(update(posi.model, .~1)$coefficients[1])*posi.model$n
  } else pull.LL <- 0

  if (class(binom.model)[1] == "speedglm") {
    # dat.all <- eval(binom.model$call$data)
    # dat.all <- cbind(dat.all, p = dat.all[,names(attr(binom.model$terms,"dataClasses"))[1]])
    if (!isTRUE(binom.model$call$fitted) || !isTRUE(posi.model$call$fitted)) {
      stop("Set fitted=TRUE in speedglm()")
    }
    dat.temp <- eval(posi.model$call$data)
    dat.temp$p <- 1
    dat.all <- merge(eval(binom.model$call$data), dat.temp, all=T)
    dat.all$p <- ifelse (dat.all$p==1,1,0)
    } else {
    dat.all <- data.frame(binom.model$data, p = binom.model$y)
  }
  Year <- dat.all[, factor[1]]
  Area <- rep(1, nrow(dat.all))
  Lon <- dat.all[,lon.name]
  Lat <- dat.all[,lat.name]
  
  obj.f <- function(x, i, k, LonLat = "Lon", criteria = criteria){
    Area <- as.numeric(Area)
    if(LonLat == "Lon"){
      # for(ii in 1:nrow(dat.all)) if((Area[ii] == k) && (Lon[ii] > x)) Area[ii] <- i+1
      Area <- sapply(1:nrow(dat.all), function(ii) {
        if ((Area[ii] == k) && (Lon[ii] > x)) i+1 else i
      })  
    }
    if(LonLat == "Lat"){
      # for(ii in 1:nrow(dat.all)) if((Area[ii] == k) && (Lat[ii] > x)) Area[ii] <- i+1
      Area <- sapply(1:nrow(dat.all), function(ii) {
        if ((Area[ii] == k) && (Lat[ii] > x)) i+1 else i
      })
    }
    
    dat.all[, factor[2]] <- as.factor(Area)
    dat.posi <- subset(dat.all, p==1)
    
    if (class(binom.model)[1]=="speedglm") { #if using "speedglm"
      res1 <- try(update(binom.model,data = dat.all,add=FALSE), silent=TRUE)
    } else {
      res1 <- try(update(binom.model,data = dat.all), silent=TRUE)
    }
    
    if (class(posi.model)[1]=="speedglm") {
      res2 <- try(update(posi.model, data = dat.posi, add=FALSE), silent=TRUE)
    } else {
      res2 <- try(update(posi.model, data = dat.posi), silent=TRUE)
    }

    if(class(res1) == "try-error" || class(res2) == "try-error") {
      obj <- Inf
    } else {
      res <- update.delta.glm(res1, res2, trend=FALSE, remove=remove, data = dat.all, pull.LL = pull.LL)
      if(criteria=="AICc") obj <- res$AICc
      if(criteria=="AIC") obj <- res$AIC
      if(criteria=="BIC") obj <- res$BIC
    }
    return(obj)
  }
  
  seq.lat1 <- seq(min.lat, max.lat, seq.lon)
  seq.lat1 <- seq.lat1[-c(1,length(seq.lat1))]
  seq.lon1 <- seq(min.lon, max.lon, seq.lat)
  seq.lon1 <- seq.lon1[-c(1,length(seq.lon1))]
  
  for(i in 1:(max.narea-1)){
    
    if (parallel) {
      cl <- makeCluster(cores, type=type)
      registerDoParallel(cl)
      exports <- c("dat.all","update.delta.glm")
      res.lon <- foreach(k=1:i, .packages = export.package, .export=exports) %dopar% {
        min.lon2 <- max(na.omit(tapply(Lon, list(Year, Area), min)[,k]))
        max.lon2 <- min(na.omit(tapply(Lon, list(Year, Area), max)[,k]))
        seq.lon2 <- seq.lon1[seq.lon1 > min.lon2 & seq.lon1 < max.lon2]
        if (min.lon2 >= max.lon2 || length(seq.lon2) == 0) res.c <- list(minimum=mean(min.lon2, max.lon2), objective=Inf) 
        else {
          obj.list <-  sapply(1:length(seq.lon2), function(ii) obj.f(x=seq.lon2[ii], i, k, LonLat="Lon", criteria=criteria))
          obj.pos <- which.min(obj.list)
          res.c <- list(minimum=seq.lon2[obj.pos], objective=obj.list[obj.pos])
        }
        res.c
      }
      stopCluster(cl)
      
      cl <- makeCluster(cores, type=type)
      registerDoParallel(cl)
      res.lat <- foreach(k =1:i, .packages = export.package, .export=exports) %dopar% {
        min.lat2 <- max(na.omit(tapply(Lat, list(Year, Area), min)[,k]))
        max.lat2 <- min(na.omit(tapply(Lat, list(Year, Area), max)[,k]))
        seq.lat2 <- seq.lat1[seq.lat1 > min.lat2 & seq.lat1 < max.lat2]
        if(min.lat2 >= max.lat2 || length(seq.lat2) == 0) res.c <- list(minimum=mean(min.lat2, max.lat2), objective=Inf) 
        else {
          obj.list <-  sapply(1:length(seq.lat2), function(ii) obj.f(x=seq.lat2[ii], i, k, LonLat="Lat", criteria=criteria))
          obj.pos <- which.min(obj.list)
          res.c <- list(minimum=seq.lat2[obj.pos], objective=obj.list[obj.pos])
        }
        res.c
      }
      stopCluster(cl)
      
    } else {
      res.lon <- lapply(1:i, function(k){
        min.lon2 <- max(na.omit(tapply(Lon, list(Year, Area), min)[,k]))
        max.lon2 <- min(na.omit(tapply(Lon, list(Year, Area), max)[,k]))
        seq.lon2 <- seq.lon1[seq.lon1 > min.lon2 & seq.lon1 < max.lon2]
        if (min.lon2 >= max.lon2 || length(seq.lon2) == 0) res.c <- list(minimum=mean(min.lon2, max.lon2), objective=Inf) 
        else {
          obj.list <-  sapply(1:length(seq.lon2), function(ii) obj.f(x=seq.lon2[ii], i, k, LonLat="Lon", criteria=criteria))
          obj.pos <- which.min(obj.list)
          res.c <- list(minimum=seq.lon2[obj.pos], objective=obj.list[obj.pos])
        }
        res.c
      })
      
      res.lat <- lapply(1:i, function(k){
        min.lat2 <- max(na.omit(tapply(Lat, list(Year, Area), min)[,k]))
        max.lat2 <- min(na.omit(tapply(Lat, list(Year, Area), max)[,k]))
        seq.lat2 <- seq.lat1[seq.lat1 > min.lat2 & seq.lat1 < max.lat2]
        if(min.lat2 >= max.lat2 || length(seq.lat2) == 0) res.c <- list(minimum=mean(min.lat2, max.lat2), objective=Inf) 
        else {
          obj.list <-  sapply(1:length(seq.lat2), function(ii) obj.f(x=seq.lat2[ii], i, k, LonLat="Lat", criteria=criteria))
          obj.pos <- which.min(obj.list)
          res.c <- list(minimum=seq.lat2[obj.pos], objective=obj.list[obj.pos])
        }
        res.c
      })
    }
    
    res.lon.all <- sapply(1:i, function(k)res.lon[[k]]$objective)
    group.lon <- which.min(res.lon.all)
    res.lon.min <- res.lon[[group.lon]]
    
    res.lat.all <- sapply(1:i, function(k)res.lat[[k]]$objective)
    group.lat <- which.min(res.lat.all)
    res.lat.min <- res.lat[[group.lat]]
    
    lonlat <- which.min(c(res.lon.min$objective, res.lat.min$objective))
    if(i >= 2) min.obj <-min(res$objective) 
    res$objective[i] <- min(res.lon.min$objective, res.lat.min$objective)
    
    if(i >= 2) Area <- as.numeric(res$area[[i-1]][,3]) else Area <- rep(1, nrow(dat.all))
    
    if(lonlat == 1){
      sep <- res.lon.min$minimum
      res$LonLat[i] <- "Lon"
      res$sep[i] <- sep
      res$sep.area[i] <- group.lon
      for(ii in 1:nrow(dat.all)) {
        if((Area[ii] == group.lon) && (Lon[ii] > sep)) Area[ii] <- i+1
      }
    }
    if(lonlat == 2){
      sep <- res.lat.min$minimum
      res$LonLat[i] <- "Lat"
      res$sep[i] <- sep
      res$sep.area[i] <- group.lat
      for(ii in 1:nrow(dat.all)){
        if((Area[ii] == group.lat) && (Lat[ii] > sep)) Area[ii] <- i+1
      }
    }
    
    dat.all <- cbind(dat.all, Area)    
    colnames(dat.all)[length(colnames(dat.all))] <- sprintf("Area%s",(i+1))
    
    dat.all[,factor[2]] <- as.factor(Area)
    dat.posi <- subset(dat.all, p == 1)
    
    if(res$objective[i] != Inf){
      if (class(binom.model)[1]=="speedglm") { #if using "speedglm"
        res1 <- try(update(binom.model,data = dat.all,add=FALSE), silent=TRUE)
      } else {
        res1 <- try(update(binom.model,data = dat.all), silent=TRUE)
      }
      
      if (class(posi.model)[1]=="speedglm") {
        res2 <- try(update(posi.model, data = dat.posi,add=FALSE), silent=TRUE)
      } else {
        res2 <- try(update(posi.model, data = dat.posi), silent=TRUE)
      }
      
      res.update <- update.delta.glm(res1, res2, trend, factor, numeric, seq.length, offset, area.w, lon.name, lat.name, min.lon, max.lon, min.lat, max.lat, seq.lon, seq.lat, scale, remove, data=dat.all,pull.LL=pull.LL)
      
      res$delta.glm[[i]] <- res.update
      res$area[[i]] <-   cbind("Lon" = Lon, "Lat" = Lat, "Area" = Area)
    }
    if(cat) cat(sprintf("[[%s]]objective:%10.2f\n",i,res$objective[i]))
      
    if((i >= 2) && (res$objective[i]-min.obj >= delta))break
  }
  
  best.pos <- which.min(res$objective)
  res$best <- res$delta.glm[[best.pos]]
  res$best$area <- res$area[[best.pos]]
  res$dat <- cbind(dat.all, "BestArea" = res$best$area[,"Area"])

  if (trend) {
    if (model.ave) {
      res$delta <- res$objective[-length(res$objective)]-min.obj
      res$delta2 <- ifelse(res$delta< delta, res$delta, 100)
      res$weight <- exp(-res$delta2/2)/sum(exp(-res$delta2/2))
      
      res$ave$p.trend0 <- rowSums(sapply(1:length(res$weight), function(j)res$weight[j]*res$delta.glm[[j]]$p.trend0))
      res$ave$m.trend0 <- rowSums(sapply(1:length(res$weight), function(j)res$weight[j]*res$delta.glm[[j]]$m.trend0))
      res$ave$y.trend0 <- rowSums(sapply(1:length(res$weight), function(j)res$weight[j]*res$delta.glm[[j]]$y.trend0))
      
      res$ave$p.trend <- res$ave$p.trend0
      res$ave$m.trend <- res$ave$m.trend0 / get(scale)(res$ave$m.trend0,na.rm=TRUE)
      res$ave$y.trend <- res$ave$y.trend0 / get(scale)(res$ave$y.trend0,na.rm=TRUE)
    }
  }
  return(res)
}




