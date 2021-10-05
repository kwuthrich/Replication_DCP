##############################################################################
# Paper: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# Figures motivating example
# DISCLAIMER: This software is provided "as is" without warranty of 
# any kind, expressed or implied.
##############################################################################

rm(list = ls())
set.seed(12345)

### Preliminaries

setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Conformal QR and DR/code/replication_package")

##############################################################################
# DGP: Y=X+X*U where X~Unif(0,1) and U~N(0,1)
##############################################################################

alpha.sig <- 0.1
n <- 5000000
x <- runif(n)
u <- rnorm(n)
xs <- seq(0.05,0.95,0.05)

# (Un)conditional length mean-based conformal prediction

uncond.reg <- 2*quantile(abs(x*u),1-alpha.sig)

# Conditional length distributional conformal prediction

cond.distr <- (qnorm(1-alpha.sig/2)-qnorm(alpha.sig/2))*xs

# Conditional coverage

crit <- quantile(abs(x*u),1-alpha.sig)
cond.cov.reg <- function(t,crit){
  pnorm(crit/t)-pnorm(-crit/t)
}

cond.cov.reg <- unlist(lapply(xs,cond.cov.reg,crit))

##############################################################################
# Figures
##############################################################################

graphics.off()
pdf("picture_intro_length.pdf", pointsize=14,width=6.0,height=6.0)
plot(c(0,1),c(0,3.5), type="n", ylab="Conditional length", xlab="X", main="(a) Conditional length 90% prediction interval")
lines(xs, cond.distr, lty=1,lwd=3)
lines(xs, rep(uncond.reg,length(xs)),lty=3,lwd=3)
legend("topleft", c("Distributional conformal prediction", "Mean-based conformal prediction"), lty = c(1,3), lwd = c(3,3),bty = "n")
dev.off()

graphics.off()
pdf("picture_intro_coverage.pdf", pointsize=14,width=6.0,height=6.0)
plot(c(0,1),c(0.6,1), type="n", ylab="Conditional coverage", xlab="X", main="(b) Conditional coverage 90% prediction interval")
lines(xs, rep(1-alpha.sig,length(xs)),lty=1,lwd=3)
lines(xs, cond.cov.reg,lty=3,lwd=3)
legend("bottomleft", c("Distributional conformal prediction", "Mean-based conformal prediction"), lty = c(1,3), lwd = c(3,3),bty = "n")
dev.off()



