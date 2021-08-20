####################################################################################
# Paper: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# Empirical illustration based on CPS data
# DISCLAIMER: This software is provided "as is" without warranty of 
# any kind, expressed or implied.
####################################################################################

# Preliminaries
rm(list = ls())
set.seed(12345)

# Packages
library(xtable)
library(hdm)
library(matrixcalc)

# Working directory
setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Conformal QR and DR/code")

# Functions
source("functions_015.R")

####################################################################################
# Data preparation
####################################################################################

# Load data
data(cps2012)

# Generate model matrix with all two-way interactions
Xo <- model.matrix(~-1 + (female + widowed + divorced + separated + nevermarried + hsd08 + hsd911
                          + hsg + cg + ad + mw + so + we + exp1 + exp2)^2, data = cps2012)

# Remove constant variables
Xo <- Xo[, which(apply(Xo, 2, var) != 0)]

# Generate wage based on log wage
Yo <- exp(cps2012$lnw)

####################################################################################
# Results
####################################################################################

# Miscoverage level
alpha.sig <- 0.1

# 20 experiments with 20% random holdout samples
T.ho <- floor(length(Yo)*0.2)

cov.mat <- leng.mat <- X.test.mat <- NULL

begin <- Sys.time()
for (r in 1:20){
  
  print(r)

  # Define training and holdout samples (while loop is to avoid singular design)
  sing <- TRUE
  while (sing==TRUE){
    
    ind <- sample(length(Yo),length(Yo),replace=FALSE)
    Y   <- Yo[ind]
    X   <- Xo[ind,]
    
    ind.test  <- (length(Y)-T.ho+1):length(Y)
    Y.test    <- Y[ind.test]
    X.test    <- cbind(X[ind.test,])
    
    obj.uneven <- aux.uneven(Y[-ind.test],cbind(X[-ind.test,]))
    
    X0 <- cbind(obj.uneven$X0)
    Y0 <- obj.uneven$Y0
    X1 <- cbind(obj.uneven$X1)
    Y1 <- obj.uneven$Y1
    
    sing  <- is.singular.matrix(t(X0)%*%X0)
  }
  
  # Choosing grids for QR and DR (we found that DR often works better with a coarser grid)
  taus    <- seq(0.001,0.999,length=200)
  ys      <- quantile(unique(c(Y0,Y1)),seq(0.001,0.999,length=200))
  
  # Applying the difference conformal prediction methods
  res.qr    <- dcp.qr(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
  res.dr    <- dcp.dr(Y0,X0,Y1,X1,Y.test,X.test,ys,taus,alpha.sig)
  res.cqr   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  res.reg   <- cp.reg(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  res.loc   <- cp.loc(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
  
  # Results
  X.test.mat <- rbind(X.test.mat,X.test)
  
  cov.mat.temp  <- cbind(res.qr$cov.qr,res.opt$cov.opt,res.dr$cov.dr,res.cqr$cov.o,res.cqr$cov.m,res.cqr$cov.r,res.reg$cov.reg,res.loc$cov.loc)
  leng.mat.temp <- cbind(res.qr$leng.qr,res.opt$leng.opt,res.dr$leng.dr,res.cqr$leng.o,res.cqr$leng.m,res.cqr$leng.r,res.reg$leng.reg,res.loc$leng.loc)

  cov.mat   <- rbind(cov.mat,cov.mat.temp)
  leng.mat  <- rbind(leng.mat,leng.mat.temp)

}
Sys.time()-begin

save.image(file=paste("RData/cps",format(Sys.time(),"%Y-%m-%d_%T"),".RData",sep=""))

####################################################################################
# Length
####################################################################################

res.leng <- colMeans(leng.mat,na.rm=TRUE)
print(xtable(rbind(res.leng),include.rownames=F,digits=2))

####################################################################################
# Coverage
####################################################################################

# Unconditional coverage
res.cov <- colMeans(cov.mat,na.rm=TRUE)
print(xtable(rbind(res.cov),include.rownames=F,digits=2))

# Estimate predicted coverage
pred.cov.qr       <- predict(glm(cov.mat[,1]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.opt      <- predict(glm(cov.mat[,2]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.dr       <- predict(glm(cov.mat[,3]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.o    <- predict(glm(cov.mat[,4]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.m    <- predict(glm(cov.mat[,5]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.cqr.r    <- predict(glm(cov.mat[,6]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.reg      <- predict(glm(cov.mat[,7]~X.test.mat,family=binomial(link="logit")),type="response")
pred.cov.loc      <- predict(glm(cov.mat[,8]~X.test.mat,family=binomial(link="logit")),type="response")

# Dispersion of conditional coverage

pred.cov <- cbind(pred.cov.qr,pred.cov.opt,pred.cov.dr,pred.cov.cqr.o,pred.cov.cqr.m,pred.cov.cqr.r,pred.cov.reg,pred.cov.loc)
# for comparison sqrt(colMeans((pred.cov-0.9)^2))*100
res.cov.cond <- apply(pred.cov,2,sd)
print(xtable(rbind(100*res.cov.cond),include.rownames=F,digits=2))

# Impose lower bound for histograms
lb.cov <- 0.5

pred.cov.qr     <- pred.cov.qr[pred.cov.qr>lb.cov]
pred.cov.opt    <- pred.cov.opt[pred.cov.opt>lb.cov]
pred.cov.dr     <- pred.cov.dr[pred.cov.dr>lb.cov]
pred.cov.cqr.o  <- pred.cov.cqr.o[pred.cov.cqr.o>lb.cov]
pred.cov.cqr.m  <- pred.cov.cqr.m[pred.cov.cqr.m>lb.cov]
pred.cov.cqr.r  <- pred.cov.cqr.r[pred.cov.cqr.r>lb.cov]
pred.cov.reg    <- pred.cov.reg[pred.cov.reg>lb.cov]
pred.cov.loc    <- pred.cov.loc[pred.cov.loc>lb.cov]

# Number of break points
breaks.cov   <- seq(lb.cov,1,length=100)

# Histograms
hg.cov.qr     <- hist(pred.cov.qr, breaks = breaks.cov,plot = FALSE)
hg.cov.opt    <- hist(pred.cov.opt, breaks = breaks.cov,plot = FALSE)
hg.cov.dr     <- hist(pred.cov.dr, breaks = breaks.cov, plot = FALSE)
hg.cov.cqr.o  <- hist(pred.cov.cqr.o, breaks = breaks.cov, plot = FALSE)
hg.cov.cqr.m  <- hist(pred.cov.cqr.m, breaks = breaks.cov, plot = FALSE)
hg.cov.cqr.r  <- hist(pred.cov.cqr.r, breaks = breaks.cov, plot = FALSE)
hg.cov.reg    <- hist(pred.cov.reg, breaks = breaks.cov, plot = FALSE)
hg.cov.loc    <- hist(pred.cov.loc, breaks = breaks.cov, plot = FALSE)

# Upper limit for y axis
upper.limit.cov <- 35000

# Generate figures
pdf("graphics/cps_cov_qr.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.qr, col = "darkgrey",axes=T, ylab="",ylim=c(0,upper.limit.cov),yaxt="n",xlab="Predicted conditional coverage",main="DCP-QR")
abline(v=0.9,lwd=3)
dev.off()
pdf("graphics/cps_cov_opt.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.opt, col = "darkgrey", ylab="",ylim=c(0,upper.limit.cov),yaxt="n",xlab="Predicted conditional coverage",main="DCP-QR*")
abline(v=0.9,lwd=3)
dev.off()
pdf("graphics/cps_cov_dr.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.dr, col = "darkgrey", ylab="",ylim=c(0,upper.limit.cov),yaxt="n", xlab="Predicted conditional coverage",main="DCP-DR")
abline(v=0.9,lwd=3)
dev.off()
pdf("graphics/cps_cov_cqro.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.cqr.o, col = "darkgrey", ylab="",ylim=c(0,upper.limit.cov),yaxt="n",xlab="Predicted conditional coverage",main="CQR")
abline(v=0.9,lwd=3)
dev.off()
pdf("graphics/cps_cov_cqrm.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.cqr.m, col = "darkgrey", ylab="",ylim=c(0,upper.limit.cov),yaxt="n",xlab="Predicted conditional coverage",main="CQR-m")
abline(v=0.9,lwd=3)
dev.off()
pdf("graphics/cps_cov_cqrr.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.cqr.r, col = "darkgrey", ylab="",ylim=c(0,upper.limit.cov),yaxt="n",xlab="Predicted conditional coverage",main="CQR-r")
abline(v=0.9,lwd=3)
dev.off()
pdf("graphics/cps_cov_reg.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.reg, col = "darkgrey", ylab="",ylim=c(0,upper.limit.cov),yaxt="n",xlab="Predicted conditional coverage",main="CP-OLS")
abline(v=0.9,lwd=3)
dev.off()
pdf("graphics/cps_cov_loc.pdf", pointsize=30,width=8.0,height=8.0)
par(mar = c(2 , 1, 1, 1))
plot(hg.cov.loc, col = "darkgrey", ylab="",ylim=c(0,upper.limit.cov),yaxt="n",xlab="Predicted conditional coverage",main="CP-loc")
abline(v=0.9,lwd=3)
dev.off()
