#######################################################################################################################################
### Parallel group multi-period cluster randomized trials 
### The effect of various measurement schemes and attrition
### Cross-sectional design
###
### author: Mirjam Moerbeek
### date last change: June 17, 2020
#######################################################################################################################################

f.varCRT=function(r,rho,m,k,R,weekdays,maxdays,omega0,gamma0,omega1,gamma1){

  ### measurement scheme
  scheme=weekdays
  if (R>1)
  {
    for(i in 1:(R-1))
      scheme=c(scheme,weekdays+i*7)
  }
  T=length(scheme)

  ### covariance matrix of responses
  distance=seq(0,R*7-1)
  covariance=toeplitz(r^distance)
  covariance=covariance[scheme,scheme]
  periods=as.factor(seq(1,(R*7)))
  Z=model.matrix(~0+periods)
  Z=Z[scheme,scheme]
  Z=Z[rep(1:(T),rep(m,T)),]
  V=diag(1-rho,m*T)+rho*Z%*%covariance%*%t(Z)
  Vinv=solve(V)

  ### design matrices for both experimental condition (fixed part of the model)
  treat=rep(0,times=nrow(Z))
  X0=cbind(Z,treat)
  treat=rep(1,times=nrow(Z))
  X1=cbind(Z,treat)
  
  ####################################################################################
  ### covariance covariance matrix of regression weights in the case of no attrition
  ####################################################################################
  if (omega0==0)
  {
    covmat=solve(k*t(X0)%*%Vinv%*%X0+k*t(X1)%*%Vinv%*%X1)
  }

  ####################################################################################
  ### covariance matrix of regression weights in the case of attrition
  ####################################################################################
  if (omega0>0){
    ### survival function 
    time=seq(0,1,length=maxdays)
    remain0 = (1 - omega0)^time^gamma0
    remain1 = (1 - omega1)^time^gamma1
    
    ### covariance matrix
    covmat <- 0
    prob0=0
    prob0.sum=0
    prob1=0
    prob1.sum=0
    for(j in (R*7):1) 
    {
      if(is.element(j,scheme)==TRUE)
      {
        scheme=scheme[scheme<=j]
        prob0=remain0[j]-prob0.sum
        prob0.sum=prob0.sum+prob0
        prob1=remain1[j]-prob1.sum
        prob1.sum=prob1.sum+prob1
        X0 <- X0[1:(length(scheme)*m),  ]
        X1 <- X1[1:(length(scheme)*m),  ]
        Z <- Z[1:(length(scheme)*m),  ]
        
        if(is.matrix(X0)==FALSE)
          X0=matrix(X0,nrow=1,ncol=length(X0))
        if(is.matrix(X1)==FALSE)
          X1=matrix(X1,nrow=1,ncol=length(X1))
        if(is.matrix(Z)==FALSE)
          Z=matrix(Z,nrow=1,ncol=length(Z))
        
        V <- diag(1-rho,length(scheme)*m)+rho*Z %*% covariance %*% t(Z)
        Vinv <- solve(V)
        covmat <- covmat + k * prob0 * t(X0) %*% Vinv %*% X0 + k * prob1 * t(X1) %*% Vinv %*% X1
      }
    }
    covmat <- solve(covmat)
  }
  return(covmat[T+1,T+1])  
}

