## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5, fig.path='Figs/',
                      echo=FALSE, warning=FALSE, message=FALSE)

library(timereg)
library(transform.hazards)

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
n1 <- 300
Samp1 <- rexp(n1,1)
cTimes <- runif(n1)*4
Samp1 <- pmin(Samp1,cTimes)
isNotCensored <- cTimes != Samp1
startTimes1 <- rep(0,n1)

aaMod1 <- aalen(Surv(startTimes1,Samp1,isNotCensored)~1)

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
dA_est1 <- diff(c(0,aaMod1$cum[,2]))
hazMatrix <- matrix(dA_est1,nrow=1)

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
F_survival <- function(X)matrix(-X,nrow=1,ncol=1)
F_survival_JacobianList <- list(function(X)matrix(-1,nrow=1,ncol=1))

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
param <- pluginEstimate(n1, hazMatrix, F_survival, F_survival_JacobianList,matrix(1,nrow=1,ncol=1),matrix(0,nrow=1,ncol=1))

## ---- fig.show='hold', echo=T-------------------------------------------------
times1 <- aaMod1$cum[,1]
plot(times1,param$X,type="s",xlim=c(0,4),xlab="t",ylab="",main="Survival")
lines(times1,param$X + 1.96 * sqrt(param$covariance[1,1,]),type="s",lty=2)
lines(times1,param$X - 1.96 * sqrt(param$covariance[1,1,]),type="s",lty=2)
lines(times1,exp(-times1),col=2)
legend("topright",c("estimated","exact"),lty=1,col=c(1,2),bty="n")

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
F_restrict <- function(X)matrix(c(-X[1],0,0,X[1]),nrow=2)
F_restrict_JacobianList <- list(function(X)matrix(c(-1,0,0,0),nrow=2),
                                function(X)matrix(c(0,0,1,0),nrow=2,byrow=T))

X0_restrict <- matrix(c(1,0),nrow=2)
V0_restrict <- matrix(0,nrow=2,ncol=2)

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
fimeTimes <- seq(0,4,length.out = 1e4+1)

tms <- sort(unique(c(fimeTimes,times1)))

hazMatrix <- matrix(0,nrow=2,ncol=length(tms))
hazMatrix[1,match(times1,tms)] <- dA_est1
hazMatrix[2,] <- diff(c(0,tms))

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
param <- pluginEstimate(n1, hazMatrix, F_restrict, F_restrict_JacobianList,X0_restrict,V0_restrict,isLebesgue = 2)

## ---- fig.show='hold', echo=T-------------------------------------------------
plot(tms,param$X[2,],type="s",xlim=c(0,4),ylim=c(0,1.5),xlab="t",ylab="",main="Restricted mean survival")
lines(tms,param$X[2,] + 1.96 * sqrt(param$covariance[2,2,]),type="s",lty=2)
lines(tms,param$X[2,] - 1.96 * sqrt(param$covariance[2,2,]),type="s",lty=2)
lines(tms,1 - exp(-tms),col=2)
legend("topleft",c("estimated","exact"),lty=1,col=c(1,2),bty="n")


## ---- fig.show='hold',results = "hide",echo=T---------------------------------
n2 <- 200
Samp2 <- rexp(n2,1.3)
censTimes2 <- runif(n2)*3
Samp2 <- pmin(Samp2,censTimes2)
isNotCensored2 <- censTimes2 != Samp2
startTimes2 <- rep(0,n2)

aaMod2 <- aalen(Surv(startTimes2,Samp2,isNotCensored2)~1)
dA_est2 <- diff(c(0,aaMod2$cum[,2]))

times2 <- aaMod2$cum[,1]
times <- sort(unique(c(times1,times2)))

hazMatrix <- matrix(0,nrow=2,ncol=length(times))
hazMatrix[1,match(times1,times)] <- dA_est1
hazMatrix[2,match(times2,times)] <- dA_est2

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
F_relsurv <- function(X)matrix(c(-X,X),nrow=1)
F_relsurv_JacobianList <- list(function(X)matrix(-1,nrow=1),
                               function(X)matrix(1,nrow=1))

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
param <- pluginEstimate(n1+n2, hazMatrix, F_relsurv, F_relsurv_JacobianList,matrix(1,nrow=1,ncol=1),matrix(0,nrow=1,ncol=1))

## ---- fig.show='hold', echo=T-------------------------------------------------
plot(times,param$X,type="s",xlim=c(0,4),ylim=c(-0.5,4),xlab="t",ylab="",main="Relative survival")
lines(times,param$X + 1.96 * sqrt(param$covariance[1,1,]),type="s",lty=2)
lines(times,param$X - 1.96 * sqrt(param$covariance[1,1,]),type="s",lty=2)
lines(times,exp(-times*(1/1.3 - 1)),col=2)
legend("topleft",c("estimated","exact"),lty=1,col=c(1,2),bty="n")


## ---- fig.show='hold',results = "hide",echo=T---------------------------------
F_cuminc <- function(X)matrix(c(-X[1],-X[1],X[1],0,0,X[1]),nrow=3,byrow=T)
F_cuminc_JacobianList <- list(function(X)matrix(c(-1,0,0,1,0,0,0,0,0),nrow=3,byrow=T),
                               function(X)matrix(c(-1,0,0,0,0,0,1,0,0),nrow=3,byrow=T))
X0_cuminc <- matrix(c(1,0,0),nrow=3)
v0_cuminc <- matrix(0,nrow=3,ncol=3)

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
# Generate the data
dfr <- data.frame(from=0,to=rexp(n1+n2,1),from.state = 1)

# Adding causes of death (to.state = 2 or 3) and censoring (to.state = 0)
dfr$to.state <- sample(c(0,2,3),n1+n2,replace=T,prob=c(0.2,0.2,0.6))


# Obtaining cause-specific cumulative hazard estimates and extracting the increments
aamod22 <- aalen(Surv(from,to,to.state %in% c(2)) ~1, data=dfr)
aamod33 <- aalen(Surv(from,to,to.state %in% c(3)) ~1, data=dfr)
times = sort(unique(c(aamod22$cum[,1],aamod33$cum[,1])))
hazMatrix = matrix(0,nrow=2,ncol=length(times))
mt1 = match(aamod22$cum[,1],times)
mt2 = match(aamod33$cum[,1],times)
hazMatrix[1,mt1] = diff(c(0, aamod22$cum[,2] ))
hazMatrix[2,mt2] = diff(c(0, aamod33$cum[,2] ))

# Transforming the estimates
param <- pluginEstimate(n1+n2,hazMatrix,F_cuminc,F_cuminc_JacobianList,X0_cuminc,v0_cuminc)

## ---- fig.show='hold', echo=T-------------------------------------------------
plot(times,param$X[2,],type="s",xlim=c(0,4),ylim=c(0,1),xlab="t",ylab="",main="Cumulative incidence")
lines(times,param$X[2,] + 1.96 * sqrt(param$covariance[2,2,]),type="s",lty=2)
lines(times,param$X[2,] - 1.96 * sqrt(param$covariance[2,2,]),type="s",lty=2)
lines(times,param$X[3,],type="s",col=3)
lines(times,param$X[3,] + 1.96 * sqrt(param$covariance[3,3,]),type="s",lty=2,col=3)
lines(times,param$X[3,] - 1.96 * sqrt(param$covariance[3,3,]),type="s",lty=2,col=3)
legend("topleft",c(expression(C^1),expression(C^2)),lty=1,col=c(1,3),bty="n")

