###################################################################
# Paper: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# Illustration optimal procedure
# DISCLAIMER: This software is provided "as is" without warranty of 
# any kind, expressed or implied.
###################################################################

# Preliminaries

rm(list = ls())
set.seed(12345)
setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Conformal QR and DR/code/replication_package")

# Draw figure

x1 <- seq(-3,3,length=100)
x2 <- seq(0,15,length=100)

y1 <- dnorm(x1,0,1)

df.chi2 <- 5
y2      <- dchisq(x2,df=df.chi2)

crit.fun <- function(tau) qchisq(tau+0.9,df=df.chi2)-qchisq(tau,df=df.chi2)

tau.ast <- optimize(crit.fun, c(0,0.1),maximum=FALSE)$minimum

round(tau.ast,digits=3) 
round(tau.ast,digits=3) + 0.9

pdf("illu_optimal_chi2.pdf", pointsize=14,width=6.0,height=6.0)
plot(x2, y2,type="l",ylab="",xlab="",yaxt="n",xaxt="n",main=expression(paste(chi^2,(5))))
polygon(c(x2[x2>=qchisq(tau.ast,df=df.chi2)&x2<=qchisq(tau.ast+0.9,df=df.chi2)],qchisq(tau.ast+0.9,df=df.chi2),qchisq(tau.ast,df=df.chi2)), c(y2[x2>=qchisq(tau.ast,df=df.chi2)&x2<=qchisq(tau.ast+0.9,df=df.chi2)],0,0), col="gray")
axis(1, at=c(qchisq(tau.ast,df=df.chi2),qchisq(tau.ast+0.9,df=df.chi2)), labels=c(expression(Q[Y](0.007)),expression(Q[Y](0.907))))
dev.off()

graphics.off()
pdf("illu_optimal_normal.pdf", pointsize=14,width=6.0,height=6.0)
plot(x1, y1, type="l",ylab="",xlab="",yaxt="n",xaxt="n",main=expression(N(0,1)))
polygon(c(x1[x1>=qnorm(0.05)&x1<=qnorm(0.95)],qnorm(0.95),qnorm(0.05)), c(y1[x1>=qnorm(0.05)&x1<=qnorm(0.95)],0,0), col="gray")
axis(1, at=c(qnorm(0.05),qnorm(0.95)), labels=c(expression(Q[Y](0.05)),expression(Q[Y](0.95))))
dev.off()


