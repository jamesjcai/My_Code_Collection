wSumTestStat <- function(Y,X){
   mU <- apply(X[Y==0,],2,sum)
   nU <- sum(Y==0)
   n <- length(Y)
   q <- (mU+1)/(2*(nU+1))
   w <- sqrt(n*q*(1-q))


   gamma<-c(rep(0,n))
   for(j in 1:n)
     gamma[j] <- sum(X[j,]/w)

   score <- sum(rank(gamma,ties.method= "min")[Y==1])

   return(score)

}

wSumTest <- function(Y,X, B=500){

    wSumT <- wSumTestStat(Y,X)
    n <- length(Y)
    wSumTs <- c(rep(0,B))
    for(i in 1:B){
       Yb <- sample(Y, n)
       wSumTs[i] <- wSumTestStat(Yb,X)
    }

    mu <- mean(wSumTs)
    sigma <- sqrt(var(wSumTs))

    if (sigma<1e-20) sigma<-1e-20
    z <- (wSumT-mu)/sigma

    c(z, 1-pchisq(z^2,df=1) )
}

############Li and Leal's CMC test:
CMC<-function(Y,X,MAF0=0.01){
    nn = length(Y)
    data1<-cbind(Y, X)
    fN<-apply(X, 2, sum)/(2*nrow(X))
    k <- sum(fN>MAF0)+1

    data3<-data1[,-1]
    data2 <- matrix(c(rep(0,k*nn)),nrow=nn)
    m=2
    for(i in 1:ncol(data3)){
       if(fN[i]>MAF0) {data2[,m] = data3[,i]; m=m+1}
       for(j in 1:nn){
           if(fN[i]<MAF0 & data3[j,i]>0) data2[j,1] = 1
       }
    }

##########combined new genotype matrix:
    Xc<-data2

nSNP=ncol(Xc)
nSample=nrow(Xc)
nCase=sum(Y==1); nControl=sum(Y==0)

XX<-as.matrix(Xc[Y==0,])
YY<-as.matrix(Xc[Y==1,])
Xbar<-apply(XX, 2, mean)
Ybar<-apply(YY, 2, mean)

Xgb<-XX
for(i in 1:nrow(XX))
    Xgb[i,]<-XX[i,]-Xbar

Ygb<-YY
for(i in 1:nrow(YY))
    Ygb[i,]<-YY[i,]-Ybar

CovS<-(t(Xgb)%*%Xgb+t(Ygb)%*%Ygb)/(nSample-2)
##gInv of S:
if (nrow(CovS)>1){
  CovS.edecomp<-eigen(CovS)
  CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
  inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
  P<-solve(CovS.edecomp$vectors)
  gInv.CovS<-t(P) %*% diag(inveigen) %*% P
  } else{ 
         if (CovS<1e-10) CovS<-1e-10; gInv.CovS<-1/CovS }

MMT<-t(Xbar-Ybar)%*%gInv.CovS%*%(Xbar-Ybar)*nCase*nControl/(nCase+nControl)
pMMT<-1-pchisq(MMT,df=nSNP)

c(MMT, pMMT)
}



# Output is the p-values of the seven tests in the order of:                 
# 1. Score 
# 2. SSU
# 3. SSUw
# 4. UminP
# 5. Sum
# 6. aSum-P
# 7. aSum

PowerUniv <- function(U,V){
      n <- dim(V)[1]

      x <- as.numeric(max(abs(U)))
      TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))
      
      return(TER)
}

UscoTest<-function(Y,X){
U<-sum((X-mean(X))*(Y-mean(Y)))
V<- mean(Y)*(1-mean(Y))*sum( (X-mean(X))^2)
if (abs(U)<1e-20) Ts<-0 else Ts<-U^2/V
if (is.na(Ts) || is.infinite(Ts) || is.nan(Ts)) Ts<-0
c(U, 1-pchisq(Ts, df=1))
}

SumTest <- function(Y,X,alpha0){
      pv<-NULL
      Uus<-NULL
      for(i in 1:ncol(X)){
         #fit  <- glm(Y~X[,i], family=binomial(logit))
         ###"if ()" added by WP, 10/22/10, since some X[,j] can be ALL 0:
         #if (is.na(fit$coefficients[-1])){
         #   beta <- c(beta, 0); pv <- c(pv,1);
         #   } else{
         #beta <- c(beta,fit$coefficients[-1])
         #pv <- c(pv,as.numeric(summary(fit)$coefficients[,4][-1]))
         #}
         ####Score test to speed up!!!
         scoRes<-UscoTest(Y,X[,i])
         Uus<-c(Uus, scoRes[1])
         pv <- c(pv, scoRes[2])
      }

      Xg <- X

      Xbar<-mean(Xg) 
      Xgb<-Xg
       
      Xgb<-Xgb-Xbar
   
      U<-t(Xg) %*% (Y-mean(Y))
      CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)

      a<-rep(1, length(U))
      a[Uus<0 & pv<alpha0] <- -1
      
      u <- sum(t(a)%*%U)
      v <- as.numeric(t(a) %*% CovS %*% (a))

      Tsum<- sum(t(a)%*%U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
      pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )       

      return(cbind(u,pTsum,v,a*U,a))
}

rSumTest <- function(Y,X,k,alpha0){
    nCase <- sum(Y==1)
    nControl <- sum(Y==0)
    n <- nCase+nControl

    u0<-pv<-v0<-u00<-NULL    
    for(sim in 1:k){
       set.seed(sim)
       case <- sample(1:n,nCase,replace = FALSE)
       Y[case] <- 1
       Y[-case] <- 0   
       fit <- SumTest(Y,X,alpha0)
       u0 <- c(u0,fit[1,1])
       pv <- c(pv,as.numeric(fit[1,2]))
       v0 <- c(v0,fit[1,3])
       u00 <- cbind(u00,fit[,4])
    }   
     
    a <- rep(1,nrow(u00))  

    mu <- mean(u0)
    v <- var(u0)

    x<-(u0-mu)^2/v
    a <- sqrt(var(x)/2)
    b <- mean(x)-a
    
    return(cbind(rep(mean(u0),k),rep(v,k),rep(a,k),rep(b,k),pv))
}


############################################
############################################
############################################



aSumTest<-function(Y, X, B=500, alpha0=0.1){

   Xg <- X

############SSU#####################
   Xbar<-apply(Xg, 2, mean) 
   Xgb<-Xg
   for(i in 1:nrow(Xg))
      Xgb[i,]<-Xg[i,]-Xbar
   
   #########SumSqUs:
   U<-t(Xg) %*% (Y-mean(Y))
   
   #SumSqU:
     Tg1<- t(U) %*% U
     ##cov of the score stats:
       CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
     ##distr of Tg1 is sum of cr Chisq_1:
       cr<-eigen(CovS, only.values=TRUE)$values
     ##approximate the distri by alpha Chisq_d + beta:
       alpha1<-sum(cr*cr*cr)/sum(cr*cr)
     beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
     d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
     alpha1<-as.real(alpha1)
     beta1<-as.real(beta1)
     d1<-as.real(d1)
     pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))

   #SumSqUw:
     diagCovS<-diag(CovS)
     diagCovS<-ifelse(diagCovS>1e-10, diagCovS, 1e-10)
     Tg2<- t(U) %*%  diag(1/diagCovS) %*% U
     ##distr of Tg1 is sum of cr Chisq_1:
       cr<-eigen(CovS %*% diag(1/diagCovS), only.values=TRUE)$values
     ##approximate the distri by alpha Chisq_d + beta:
       alpha2<-sum(cr*cr*cr)/sum(cr*cr)
       beta2<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
       d2<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
       alpha2<-as.real(alpha2)
       beta2<-as.real(beta2)
       d2<-as.real(d2)
     pTg2<-as.numeric(1-pchisq((Tg2-beta2)/alpha2, d2))


  ##########score test:
    ##gInv of CovS:
      CovS.edecomp<-eigen(CovS)
      CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
      inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
      P<-solve(CovS.edecomp$vectors)
      gInv.CovS<-t(P) %*% diag(inveigen) %*% P
    Tscore<- t(U) %*% gInv.CovS  %*% U
    pTscore<-as.numeric( 1-pchisq(Tscore, CovS.rank) )

  #### V1 %*% t(V1)=CovS:
  #  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  #  V1<-CovS.edecomp$vectors %*% diag(sqrt(CovS.ev))

  ##univariate/marginal tests:
    
    Tus<-as.vector(abs(U)/sqrt(diagCovS))

    ####added by WP: 10/24/10:
    nSNP=length(U)

    Vs <- matrix(c(rep(0,nSNP^2)),nrow=nSNP)
    for(i in 1:nSNP){
      for(j in 1:nSNP){
          if (abs(CovS[i,j])>1e-20)
             Vs[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
          else   Vs[i,j] <- 1e-20
        }
      }
    pTus <- as.numeric(PowerUniv(Tus,Vs))



##########SumTest########################
      a<-rep(1, length(U))
      Tsum<- sum(t(a)%*%U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
      pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )  


##########aSumTest
     fit0 <- SumTest(Y,Xg,alpha0)
     u <- fit0[1,1]
     pv <- fit0[1,2]
     v <- fit0[1,3]
     fit <- rSumTest(Y,Xg,B,alpha0)
     u0 <- fit[1,1] 
     v0 <- fit[1,2]
     a <- fit[1,3]
     b <- fit[1,4]
     pv0 <- fit[,5]

   aSumP <- sum(pv>pv0)/length(pv0)
   aSum <- as.numeric( 1-pchisq(abs(((u-u0)^2/v0-b)/a), 1) )
   
   return(c(pTscore, pTg1, pTg2, pTus, pTsum, aSumP, aSum))
}  








