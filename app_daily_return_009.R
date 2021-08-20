###################################################################################
# Paper: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# Empirical illustration: predicting daily returns
# DISCLAIMER: This software is provided "as is" without warranty of 
# any kind, expressed or implied.
###################################################################################

# Preliminaries
rm(list = ls())
set.seed(12345)

# Working directory
setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Conformal QR and DR/code")
#setwd("D:/Dropbox/Co-authors/SC LASSO with Kaspar and Victor/Conformal QR and DR/code")

# Functions
source("functions_015.R")

###################################################################################
# Data 
# Source: Kenneth R. French data library
# Accessed: August 17, 2021
# URL: http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
###################################################################################

# Load data
data.daily <- read.csv("data_us_08172021.csv")

# Generate lagged realized volatility
data.daily$Mkt <- data.daily$MktRF+data.daily$RF
Tt <- nrow(data.daily)
data.daily$lrealvol <- NA 
for (t in 23:Tt){
  data.daily$lrealvol[t] <- sqrt(sum(data.daily$Mkt[(t-22):(t-1)]^2))
}

# Omit lags from computing realized volatility
data.analysis <- data.daily[23:Tt,]

Yo <- data.analysis$Mkt
Xo <- as.matrix(data.analysis$lrealvol)

To <- length(Yo)

###################################################################################
# Analysis
###################################################################################

# Miscoverage level
alpha.sig <- 0.1 

# 5 consecutive exercises with a 10% holdout sample for testing
T.ho <- floor(length(Yo)*(0.10))

cov.mat <- leng.mat <- X.test.vec <-  NULL

begin <- Sys.time()
for (r in 1:5){
  
  # Define training and holdout samples
  ind.cp    <- ((r-1)*T.ho+1):(floor(To*(0.50))+(r-1)*T.ho+1)
  ind.test  <- (max(ind.cp)+1):(max(ind.cp)+T.ho)

  print(c(min(ind.cp),max(ind.cp)))
  print(c(min(ind.test),max(ind.test)))
  
  Y.test  <- Yo[ind.test]
  X.test  <- cbind(Xo[ind.test,])

  obj.uneven <- aux.uneven(Yo[ind.cp],cbind(Xo[ind.cp,]))
  
  X0 <- cbind(obj.uneven$X0)
  Y0 <- obj.uneven$Y0
  X1 <- cbind(obj.uneven$X1)
  Y1 <- obj.uneven$Y1

  # Applying the difference conformal prediction methods
  taus    <- seq(0.001,0.999,length=200)
  ys      <- quantile(unique(c(Y0,Y1)),seq(0.001,0.999,length=200))
  
  # Applying the difference conformal prediction methods
  res.qr    <- dcp.qr(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.dr    <- dcp.dr(Y0,X0,Y1,X1,Y.test,X.test,ys,taus,alpha.sig)
  res.cqr   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  res.reg   <- cp.reg(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  res.loc   <- cp.loc(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  
  X.test.vec <- c(X.test.vec,X.test)
  
  # Return results
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.dr$cov.dr,res.cqr$cov.o,
                         res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg,res.loc$cov.loc)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.dr$leng.dr,res.cqr$leng.o,
                         res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg,res.loc$leng.loc)
  
  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)
  
}
Sys.time()-begin

save.image(file=paste("RData/return",format(Sys.time(),"%Y-%m-%d_%T"),".RData",sep=""))

# Average length
colMeans(cov.mat)
colMeans(leng.mat)

# Conditional length and conditional coverage
num.seg <- 20
cov.cond  <- binning(X.test.vec,cov.mat,num.seg)$cond
leng.cond <- binning(X.test.vec,leng.mat,num.seg)$cond

###################################################################################
# Figures
###################################################################################

graphics.off()
pdf("graphics/conditional_coverage_returns.pdf",pointsize=22,width=16.0,height=9.0)
plot(c(1,num.seg),c(0.4,1), type="n", ylab="Conditional coverage", xlab="Bin", main="Conditional coverage 90% prediction intervals")

lines((1:num.seg),cov.cond[,4],type="b",pch=4,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),cov.cond[,5],type="b",pch=5,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),cov.cond[,6],type="b",pch=6,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),cov.cond[,7],type="b",pch=8,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),cov.cond[,8],type="b",pch=11,lwd=3,col=rgb(1,0.8,0))

lines((1:num.seg),cov.cond[,1],type="b",pch=1,lwd=3,col=rgb(0,0,0.8))
lines((1:num.seg),cov.cond[,2],type="b",pch=2,lwd=3,col=rgb(0,0,0.8))
lines((1:num.seg),cov.cond[,3],type="b",pch=3,lwd=3,col=rgb(0,0,0.8))

abline(h=0.9,lty=3,lwd=1)
legend("bottomleft", c("DCP-QR","DCP-QR*", "DCP-DR","CQR","CQR-m","CQR-r","CP-OLS","CP-loc"), pch=c(1,2,3,4,5,6,8,11), 
       lwd = c(4,4,4,4,4,4,4,4),col=c(rgb(0,0,0.8),rgb(0,0,0.8),rgb(0,0,0.8),rgb(1,0.8,0),rgb(1,0.8,0), rgb(1,0.8,0),rgb(1,0.8,0),rgb(1,0.8,0)),bty = "n")
dev.off()

graphics.off()
pdf("graphics/conditional_length_returns.pdf",pointsize=22,width=16.0,height=9.0)
plot(c(1,num.seg),c(0,10), type="n", ylab="Conditional length", xlab="Bin", main="Conditional length 90% prediction intervals")

lines((1:num.seg),leng.cond[,4],type="b",pch=4,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),leng.cond[,5],type="b",pch=5,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),leng.cond[,6],type="b",pch=6,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),leng.cond[,7],type="b",pch=8,lwd=3,col=rgb(1,0.8,0))
lines((1:num.seg),leng.cond[,8],type="b",pch=11,lwd=3,col=rgb(1,0.8,0))

lines((1:num.seg),leng.cond[,1],type="b",pch=1,lwd=3,col=rgb(0,0,0.8))
lines((1:num.seg),leng.cond[,2],type="b",pch=2,lwd=3,col=rgb(0,0,0.8))
lines((1:num.seg),leng.cond[,3],type="b",pch=3,lwd=3,col=rgb(0,0,0.8))

legend("topleft", c("DCP-QR","DCP-QR*", "DCP-DR","CQR","CQR-m","CQR-r","CP-OLS","CP-loc"), pch=c(1,2,3,4,5,6,8,11), 
       lwd = c(4,4,4,4,4,4,4,4),col=c(rgb(0,0,0.8),rgb(0,0,0.8),rgb(0,0,0.8),rgb(1,0.8,0),rgb(1,0.8,0), rgb(1,0.8,0),rgb(1,0.8,0),rgb(1,0.8,0)),bty = "n")
dev.off()

