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
legend("topright",c("estimated","true"),lty=1,col=c(1,2),bty="n")

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
F_restrict <- function(X)matrix(c(-X[1],0,0,X[1]),nrow=2)
F_restrict_JacobianList <- list(function(X)matrix(c(-1,0,0,0),nrow=2),
                                function(X)matrix(c(0,0,1,0),nrow=2,byrow=T))

X0_restrict <- matrix(c(1,0),nrow=2)
V0_restrict <- matrix(0,nrow=2,ncol=2)

## ---- fig.show='hold',results = "hide",echo=T---------------------------------
fineTimes <- seq(0,4,length.out = 1e4+1)

tms <- sort(unique(c(fineTimes,times1)))

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
legend("topleft",c("estimated","true"),lty=1,col=c(1,2),bty="n")


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
legend("topleft",c("estimated","true"),lty=1,col=c(1,2),bty="n")


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
tmss = sort(unique(c(aamod22$cum[,1],aamod33$cum[,1])))
hazMatrix = matrix(0,nrow=2,ncol=length(tmss))
mt1 = match(aamod22$cum[,1],tmss)
mt2 = match(aamod33$cum[,1],tmss)
hazMatrix[1,mt1] = diff(c(0, aamod22$cum[,2] ))
hazMatrix[2,mt2] = diff(c(0, aamod33$cum[,2] ))

# Transforming the estimates
param <- pluginEstimate(n1+n2,hazMatrix,F_cuminc,F_cuminc_JacobianList,X0_cuminc,v0_cuminc)

## ---- fig.show='hold', echo=T-------------------------------------------------
plot(tmss,param$X[2,],type="s",xlim=c(0,4),ylim=c(0,1),xlab="t",ylab="",main="Cumulative incidence")
lines(tmss,param$X[2,] + 1.96 * sqrt(param$covariance[2,2,]),type="s",lty=2)
lines(tmss,param$X[2,] - 1.96 * sqrt(param$covariance[2,2,]),type="s",lty=2)
lines(tmss,param$X[3,],type="s",col=3)
lines(tmss,param$X[3,] + 1.96 * sqrt(param$covariance[3,3,]),type="s",lty=2,col=3)
lines(tmss,param$X[3,] - 1.96 * sqrt(param$covariance[3,3,]),type="s",lty=2,col=3)
legend("topleft",c(expression(C^1),expression(C^2)),lty=1,col=c(1,3),bty="n")

## ---- fig.show='hold', echo=T-------------------------------------------------
F_restrict_diff <- function(X)matrix(c(X[4],0,0,
                                       X[5],0,0,
                                       X[4]-X[5],0,0,
                                       0,-X[4],0,
                                       0,0,-X[5]),nrow=5,byrow=T)
F_restrict_diff_JacobianList <- list(function(X)matrix(c(0,0,0,1,0,
                                                         0,0,0,0,1,
                                                         0,0,0,1,-1,
                                                         rep(0,10)),nrow=5,byrow=T),
                                     function(X)matrix(c(rep(0,18),-1,rep(0,6)),nrow=5,byrow=T),
                                     function(X)matrix(c(rep(0,24),-1),nrow=5,byrow=T))

X0_restrict_diff <- matrix(c(0,0,0,1,1),nrow=5)
V0_restrict_diff <- matrix(0,nrow=5,ncol=5)

## ---- fig.show='hold', echo=T-------------------------------------------------

aa_male = aalen(Surv(futime,death==1)~1, data = flchain[flchain$sex=="M",]) 
aa_female = aalen(Surv(futime,death==1)~1, data = flchain[flchain$sex=="F",]) 

tms <- sort(unique(c(fineTimes,aa_male$cum[,1], aa_female$cum[,1])))
mt_male = match(aa_male$cum[,1],tms)
mt_female = match(aa_female$cum[,1],tms)

hazMatrix <- matrix(0,nrow=3,ncol=length(tms))
hazMatrix[1,] <- diff(c(0,tms))
hazMatrix[2,mt_male] <- diff(c(0,aa_male$cum[,2]))
hazMatrix[3,mt_female] <- diff(c(0,aa_female$cum[,2]))
  
  
plugEst <- pluginEstimate(1, hazMatrix, F_restrict_diff, F_restrict_diff_JacobianList,X0_restrict_diff,V0_restrict_diff,isLebesgue = 1)
param = list(plugEst=plugEst, tms=tms)

## ---- fig.show='hold', echo=T-------------------------------------------------
ylims = c(1.1*min(param$plugEst$X[3,] - 1.96 * sqrt(param$plugEst$covariance[3,3,])),
          1.1*max(param$plugEst$X[3,] + 1.96 * sqrt(param$plugEst$covariance[3,3,])))
plot(param$tms,param$plugEst$X[3,],type="s",ylim=ylims,ylab="",main="RMS difference",xlab="years")
lines(param$tms,param$plugEst$X[3,] + 1.96 * sqrt(param$plugEst$covariance[3,3,]),type="s",lty=2)
lines(param$tms,param$plugEst$X[3,] - 1.96 * sqrt(param$plugEst$covariance[3,3,]),type="s",lty=2)
abline(a=0,b=0,col=2)

## ---- fig.show='hold', echo=T-------------------------------------------------
F_cuminc_diff <- function(X)matrix( c(-X[1],-X[1],0,0,
                                      0,0,-X[2],-X[2],
                                      X[1],0,0,0,
                                      0,X[1],0,0,
                                      0,0,X[2],0,
                                      0,0,0,X[2],
                                      X[1],0,-X[2],0,
                                      0,X[1],0,-X[2]) ,nrow=8,byrow=T)
F_cuminc_diff_JacobianList <- list(function(X)matrix(c(-1, rep(0,7),
                                                       rep(0,8),
                                                       1, rep(0,7),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       1, rep(0,7),
                                                       rep(0,8)),nrow=8,byrow=T),
                                   function(X)matrix(c(-1, rep(0,7),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       1, rep(0,7),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       1, rep(0,7)),nrow=8,byrow=T),
                                   function(X)matrix(c(rep(0,8),
                                                       0,-1, rep(0,6),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       0,1, rep(0,6),
                                                       rep(0,8),
                                                       0,-1, rep(0,6),
                                                       rep(0,8)),nrow=8,byrow=T),
                                   function(X)matrix(c(rep(0,8),
                                                       0,-1, rep(0,6),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       rep(0,8),
                                                       0,1, rep(0,6),
                                                       rep(0,8),
                                                       0,-1, rep(0,6)),nrow=8,byrow=T))


X0_cuminc_diff <- matrix(c(1,1,0,0,0,0,0,0),nrow=8)
v0_cuminc_diff <- matrix(0,nrow=8,ncol=8)

## ---- fig.show='hold', echo=T-------------------------------------------------

mgus2$etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
mgus2$event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))


fr_M <- mgus2[mgus2$sex == "M",]
fr_M$time <- fr_M$etime/12
fr_F <- mgus2[mgus2$sex == "F",]
fr_F$time <- fr_F$etime/12

fit_PCM_M = aalen(Surv(time,event=="pcm")~1,data=fr_M)
fit_death_M = aalen(Surv(time,event=="death")~1,data=fr_M)

fit_PCM_F = aalen(Surv(time,event=="pcm")~1,data=fr_F)
fit_death_F = aalen(Surv(time,event=="death")~1,data=fr_F)

tms2 = sort(unique(c(fit_PCM_M$cum[,1],fit_PCM_F$cum[,1],
                     fit_death_M$cum[,1],fit_death_F$cum[,1])))

dA_PCM_M = dA_death_M = dA_PCM_F = dA_death_F = rep(0, length(tms2))

dA_PCM_M[match(fit_PCM_M$cum[,1], tms2)] = diff(c(0,fit_PCM_M$cum[,2]))
dA_PCM_F[match(fit_PCM_F$cum[,1], tms2)] = diff(c(0,fit_PCM_F$cum[,2]))

dA_death_M[match(fit_death_M$cum[,1], tms2)] = diff(c(0,fit_death_M$cum[,2]))
dA_death_F[match(fit_death_F$cum[,1], tms2)] = diff(c(0,fit_death_F$cum[,2]))

hazMatrix <- matrix(0,nrow=4,ncol=length(tms2))
hazMatrix[1,] <- dA_PCM_M
hazMatrix[2,] <- dA_death_M
hazMatrix[3,] <- dA_PCM_F
hazMatrix[4,] <- dA_death_F

plugEst <- pluginEstimate(1, hazMatrix, F_cuminc_diff, F_cuminc_diff_JacobianList,X0_cuminc_diff,v0_cuminc_diff)
param = list(plugEst=plugEst, tms=tms2)

## ---- fig.show='hold', echo=T-------------------------------------------------

# Colors
male_color <- "red"
female_color <- "blue"
diff_color <- "gray"

par(mfrow=c(1,2))
# Plot with colors and labels
plot(param$tms, param$plugEst$X[3,],type="s",ylim=c(-0.2,0.9),main="Cum.inc. PCM",
     ylab="",xlab="Years", col=male_color)
lines(param$tms, param$plugEst$X[3,] + 1.96 * sqrt(param$plugEst$covariance[3,3,]),type="s",lty=2, col=male_color)
lines(param$tms, param$plugEst$X[3,] - 1.96 * sqrt(param$plugEst$covariance[3,3,]),type="s",lty=2, col=male_color)

lines(param$tms, param$plugEst$X[5,],type="s", col=female_color)
lines(param$tms, param$plugEst$X[5,] + 1.96 * sqrt(param$plugEst$covariance[5,5,]),type="s",lty=2, col=female_color)
lines(param$tms, param$plugEst$X[5,] - 1.96 * sqrt(param$plugEst$covariance[5,5,]),type="s",lty=2, col=female_color)

lines(param$tms, param$plugEst$X[7,],type="s", col=diff_color)
lines(param$tms, param$plugEst$X[7,] + 1.96 * sqrt(param$plugEst$covariance[7,7,]),type="s",lty=2, col=diff_color)
lines(param$tms, param$plugEst$X[7,] - 1.96 * sqrt(param$plugEst$covariance[7,7,]),type="s",lty=2, col=diff_color)

abline(a=0,b=0,col=2)

# Add legend
legend("topleft", legend = c("Males", "Females", "Difference"), col = c(male_color, female_color, diff_color), lty = 1,bty="n")




# Plot with colors and labels
plot(param$tms, param$plugEst$X[4,],type="s",ylim=c(-0.2,0.9),main="Cum.inc. death",
     ylab="",xlab="Years", col=male_color)
lines(param$tms, param$plugEst$X[4,] + 1.96 * sqrt(param$plugEst$covariance[4,4,]),type="s",lty=2, col=male_color)
lines(param$tms, param$plugEst$X[4,] - 1.96 * sqrt(param$plugEst$covariance[4,4,]),type="s",lty=2, col=male_color)

lines(param$tms, param$plugEst$X[6,],type="s", col=female_color)
lines(param$tms, param$plugEst$X[6,] + 1.96 * sqrt(param$plugEst$covariance[6,6,]),type="s",lty=2, col=female_color)
lines(param$tms, param$plugEst$X[6,] - 1.96 * sqrt(param$plugEst$covariance[6,6,]),type="s",lty=2, col=female_color)

lines(param$tms, param$plugEst$X[8,],type="s", col=diff_color)
lines(param$tms, param$plugEst$X[8,] + 1.96 * sqrt(param$plugEst$covariance[8,8,]),type="s",lty=2, col=diff_color)
lines(param$tms, param$plugEst$X[8,] - 1.96 * sqrt(param$plugEst$covariance[8,8,]),type="s",lty=2, col=diff_color)

abline(a=0,b=0,col=2)

# Add legend
legend("topleft", legend = c("Males", "Females", "Difference"), col = c(male_color, female_color, diff_color), lty = 1,bty="n")



