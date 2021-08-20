#############################################################################
# Functions: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu 
#############################################################################

library(quantreg)

### Auxiliary functions

aux.uneven <- function(Y,X){
  T01     <- length(Y)
  T01even <- floor(T01/2)*2
  
  Yeven <- Y[(T01-T01even+1):T01]
  Xeven <- cbind(X[(length(Y)-T01+1):T01,])
  
  ind0 <- 1:(T01even/2)
  ind1 <- (T01even/2+1):T01even
  
  Y0 <- Yeven[ind0]
  X0 <- cbind(Xeven[ind0,])
  Y1 <- Yeven[ind1]
  X1 <- cbind(Xeven[ind1,])
  
  return(list(Y0=Y0,X0=X0,Y1=Y1,X1=X1))
  
}

### Methods

# DCP-QR
dcp.qr <- function(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig){
  
  beta.qr <- matrix(NA,dim(X0)[2]+1,length(taus))
  for (t in 1:length(taus)){
    beta.qr[,t] <- rq.fit.br(cbind(1,X0),Y0,tau=taus[t])$coefficients
  }
  tQ.yx   <- cbind(1,X1)%*%beta.qr
  Q.yx    <- t(apply(tQ.yx,1,FUN=sort))
  u.hat   <- rowMeans((Q.yx <= matrix(Y1,length(Y1),length(taus))))
  cs      <- abs(u.hat-0.5)

  tQ.test   <- cbind(1,X.test)%*%beta.qr
  Q.test    <- t(apply(tQ.test,1,FUN=sort))
  u.test    <- rowMeans((Q.test  <= matrix(Y.test,length(Y.test),length(taus))))
  cs.test   <- abs(u.test-0.5)
  
  k         <- ceiling((1-alpha.sig)*(1+length(Y1)))
  threshold <- sort(cs)[k]
  
  ci.grid   <- abs(taus - 0.5)
  cov.qr    <- (cs.test <= threshold)
  
  lb <- ub <- rep(NA,length(Y.test))
  for (i in 1:length(Y.test)){
    ci     <- Q.test[i,(ci.grid <= threshold)]
    ub[i]  <- max(ci)
    lb[i]  <- min(ci)
  }
  
  leng.qr <- ub-lb
  leng.qr[which(leng.qr==-Inf)] <- NA
  
  return(list(cov.qr=cov.qr,leng.qr=leng.qr))
  
}

# DCP-QR*
dcp.opt <- function(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig){

  XXX <- rbind(X1,X.test)
  YYY <- c(Y1,Y.test)
  
  beta.qr <- matrix(NA,dim(X0)[2]+1,length(taus))
  for (t in 1:length(taus)){
    beta.qr[,t] <- rq.fit.br(cbind(1,X0),Y0,tau=taus[t])$coefficients
  }
  tQ.yx   <- cbind(1,XXX)%*%beta.qr
  Q.yx    <- t(apply(tQ.yx,1,FUN=sort))
  u.hat <- rowMeans((Q.yx <= matrix(YYY,length(YYY),dim(beta.qr)[2])))
  
  bhat <- rep(NA,length(YYY))
  b.grid <- taus[taus<=alpha.sig]
  for (t in 1:length(YYY)){
    leng <- rep(NA,length(b.grid))
    leng.test <- rep(NA,length(b.grid))
    for (b in 1:length(b.grid)){
      Q.yx.u <- approx(x=taus,y=Q.yx[t,],xout=(b.grid[b]+1-alpha.sig),rule=2)$y
      leng[b] <- Q.yx.u -Q.yx[t,b]
    }
    bhat[t] <- b.grid[which.min(leng)]
  }
  
  ind.test <- (length(Y1)+1):length(YYY)
  
  cs.opt <- abs(u.hat-bhat-(1-alpha.sig)/2)
  
  k           <- ceiling((1-alpha.sig)*(1+length(Y1)))
  threshold   <- sort(cs.opt[-ind.test])[k]
  
  cov.opt   <- (cs.opt[ind.test] <= threshold)
  
  leng.opt <- NULL
  for (t in ind.test){
    ci.grid <- abs(taus - bhat[t]-(1-alpha.sig)/2)
    ci <- Q.yx[t,(ci.grid <= threshold)]
    ub <- max(ci)
    lb <- min(ci)
    leng.opt <- c(leng.opt,ub-lb)
  }
  
  leng.opt[which(leng.opt==-Inf)] <- NA
  
  return(list(cov.opt=cov.opt,leng.opt=leng.opt))  
  
}

# DCP-DR
dcp.dr <- function(Y0,X0,Y1,X1,Y.test,X.test,ys,taus,alpha.sig){
  
  beta.dr <- matrix(NA,dim(X0)[2]+1,length(ys))
  for (y in 1:length(ys)){
    beta.dr[,y] <- glm.fit(cbind(1,X0),(Y0<=ys[y]),family=binomial(link="logit"))$coefficients
  }
  tF.yx   <- plogis(cbind(1,X1)%*%beta.dr)
  F.yx    <- t(apply(tF.yx,1,FUN=sort))
  
  cs <- rep(NA,length(Y1))
  for (t in 1:length(Y1)){
    u.hat <- approx(x=ys,y=F.yx[t,],xout=Y1[t],rule=2)$y
    cs[t] <- abs(u.hat-0.5)
  }
  
  tF.test   <- plogis(cbind(1,X.test)%*%beta.dr)
  F.test    <- t(apply(tF.test,1,FUN=sort))
  
  cs.test <- rep(NA,length(Y.test))
  for (t in 1:length(Y.test)){
    u.hat <- approx(x=ys,y=F.test[t,],xout=Y.test[t],rule=2)$y
    cs.test[t] <- abs(u.hat-0.5)
  }
  
  k           <- ceiling((1-alpha.sig)*(1+length(Y1)))
  threshold   <- sort(cs)[k]
  
  ci.grid   <- abs(taus - 0.5)
  cov.dr    <- (cs.test <= threshold)
  
  lb <- ub <- rep(NA,length(Y.test))
  for (i in 1:length(Y.test)){
    ci     <- ys[abs(F.test[i,]-0.5)<=threshold]
    ub[i]  <- max(ci)
    lb[i]  <- min(ci)
  }
  
  leng.dr <- ub-lb
  
  leng.dr[which(leng.dr==-Inf)] <- NA
  
  return(list(cov.dr=cov.dr,leng.dr=leng.dr))
  
}

# CQR, CQR-m, CQR-r
cqr <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
  
  beta.lo <- rq.fit.br(cbind(1,X0),Y0,tau=(alpha.sig/2))$coefficients
  beta.hi <- rq.fit.br(cbind(1,X0),Y0,tau=(1-alpha.sig/2))$coefficients
  beta.50 <- rq.fit.br(cbind(1,X0),Y0,tau=0.5)$coefficients

  tq.lo  <- cbind(1,X1)%*%beta.lo
  tq.hi  <- cbind(1,X1)%*%beta.hi
  tq.50  <- cbind(1,X1)%*%beta.50
  
  qsr <- t(apply(cbind(tq.lo,tq.50,tq.hi),1,FUN=sort))
  
  q.lo <- qsr[,1]
  q.50 <- qsr[,2]
  q.hi <- qsr[,3]  

  Eo.vec <- Em.vec <- Er.vec <- rep(NA,length(Y1))
  for (t in 1:(length(Y1))){
    Eo.vec[t]   <- max(q.lo[t]-Y1[t],Y1[t]-q.hi[t])
    Em.vec[t]   <- max((q.lo[t]-Y1[t])/(q.50[t]-q.lo[t]),(Y1[t]-q.hi[t])/(q.hi[t]-q.50[t]))
    Er.vec[t]   <- max((q.lo[t]-Y1[t])/(q.hi[t]-q.lo[t]),(Y1[t]-q.hi[t])/(q.hi[t]-q.lo[t]))
  }
  
  k     <- ceiling((1-alpha.sig)*(1+length(Y1)))
  Q.Eo  <- sort(Eo.vec)[k]
  Q.Em  <- sort(Em.vec)[k]
  Q.Er  <- sort(Er.vec)[k]

  tq.test.lo <- cbind(1,X.test)%*%beta.lo
  tq.test.50 <- cbind(1,X.test)%*%beta.50
  tq.test.hi <- cbind(1,X.test)%*%beta.hi
  
  qs.test <- t(apply(cbind(tq.test.lo,tq.test.50,tq.test.hi),1,sort))
  
  q.test.lo <- qs.test[,1]
  q.test.50 <- qs.test[,2]
  q.test.hi <- qs.test[,3] 
  
  lb.o  <- q.test.lo - Q.Eo
  ub.o  <- q.test.hi + Q.Eo
  lb.m  <- q.test.lo - Q.Em * (q.test.50-q.test.lo) 
  ub.m  <- q.test.hi + Q.Em * (q.test.hi-q.test.50)  
  lb.r  <- q.test.lo - Q.Er * (q.test.hi-q.test.lo)  
  ub.r  <- q.test.hi + Q.Er * (q.test.hi-q.test.lo) 
  
  cov.o <- (Y.test<=ub.o & Y.test>=lb.o)
  cov.m <- (Y.test<=ub.m & Y.test>=lb.m)
  cov.r <- (Y.test<=ub.r & Y.test>=lb.r)
  
  leng.o <- ub.o-lb.o
  leng.m <- ub.m-lb.m
  leng.r <- ub.r-lb.r
  
  return(list(cov.o=cov.o,cov.m=cov.m,cov.r=cov.r,leng.o=leng.o,leng.m=leng.m,leng.r=leng.r))
  
}

# CP-OLS

cp.reg <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
  beta.reg    <- cbind(lm.fit(cbind(1,X0),Y0)$coefficients)
  cs          <- abs(Y1-cbind(1,X1)%*%beta.reg)
  k           <- ceiling((1-alpha.sig)*(1+length(Y1)))
  threshold   <- sort(cs)[k]
  cov.reg     <- (abs(Y.test-cbind(1,X.test)%*%beta.reg) <= threshold)
  leng.reg    <- rep(2*threshold,length(Y.test))
  return(list(cov.reg=cov.reg,leng.reg=leng.reg))
}

# CP-loc (we take absolute values of the MAD predictions to avoid negative values; 
# this worked better in simulations than adding a trimming constant)

cp.loc <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
  beta.reg    <- cbind(lm.fit(cbind(1,X0),Y0)$coefficients)
  absR0       <- abs(Y0-cbind(1,X0)%*%beta.reg)
  beta.sig    <- cbind(lm.fit(cbind(1,X0),absR0)$coefficients)
  sig1        <- abs(cbind(1,X1)%*%beta.sig)
  absR1       <- abs(Y1-cbind(1,X1)%*%beta.reg)
  cs          <- absR1/sig1
  k           <- ceiling((1-alpha.sig)*(1+length(Y1)))
  threshold   <- sort(cs)[k]
  lb          <- cbind(1,X.test)%*%beta.reg - threshold*abs(cbind(1,X.test)%*%beta.sig)
  ub          <- cbind(1,X.test)%*%beta.reg + threshold*abs(cbind(1,X.test)%*%beta.sig) 
  cov.loc     <- Y.test <= ub & Y.test >=lb
  leng.loc    <- ub - lb
  return(list(cov.loc=cov.loc,leng.loc=leng.loc))
}


### Graphs based on quantile bins

binning <- function(X,res.mat,num.seg){
  cond      <- matrix(NA,num.seg,dim(res.mat)[2])
  xs        <- matrix(NA,num.seg,1)
  quantiles <- seq(0,1,length=num.seg+1)
  for (i in 2:(num.seg+1)){
    q1                <- quantile(X,quantiles[i])
    q0                <- quantile(X,quantiles[i-1])
    ind               <- which((X<=q1)*(X>q0)==1)
    cond[(i-1),]      <- colMeans(res.mat[ind,],na.rm=TRUE)
    xs[(i-1),]        <- (q1+q0)/2 # midpoint
  }
  return(list(xs=xs,cond=cond))
}


