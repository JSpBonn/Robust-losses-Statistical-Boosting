# Huber:
AdaptHuber <- function(tau= 0.8213748 , k=NULL , pi=NULL , q = 0.5){ # 95% efficiency: k=1.345 translated to tau=0.8213748
  k_to_tau <- function(k){
    return((pnorm(k)-pnorm(-k))) # =1-2*pnorm(-k)
  }
  
  k_to_pi <- function(k){
    return((1/(1-2*pnorm(-k)+2*dnorm(k)/k)-1)*(-1))
  }
  if(tau!=0.8213748){
    if((0<=tau)*(tau<=1)){message(paste("using tau = ",tau,sep=""))
      tau_to_k <- function(tau){
        to_minimize2 <- function(k){ return(abs(pnorm(k)-pnorm(-k)-tau))}
        return(optimize(f=to_minimize2,interval = c(0.1,8.0),maximum=FALSE)$minimum)
      }
      k=tau_to_k(tau)
      message(paste("corresponding to k = ",k," for standard gaussian distribution, recommanded is k in [1,2] for robust results",sep=""))
      pi=k_to_pi(k)
      message(paste("corresponding to pi = ", pi,sep=""))
    } else {stop("tau is a quantile, choose value in [0,1]")}
  } else{
    if(!is.null(k)){
      if(k>0){message(paste("using k = ",k," for standard gaussian distribution, recommanded is k in [1,2] for robust results",sep=""))
        pi=k_to_pi(k)
        tau=k_to_tau(k)
        message(paste("corresponding to pi = ",pi,sep=""))
        message(paste("transferred to tau = ",tau, sep=""))
      } else {stop("k has to be a positive value")}
    } else{
      if(!is.null(pi)){
        if((0<=pi)*(pi<=1)){  message(paste("using pi = ",pi,sep=""))
          pi_to_k <- function(pi){
            to_minimize<- function(k){return(abs((1/(1-2*pnorm(-k)+2*dnorm(k)/k)-1)*(-1)-pi))}
            return(optimize(f=to_minimize,interval = c(0.01,8.00),maximum=FALSE)$minimum)
          }
          k=pi_to_k(pi)
          tau=k_to_tau(k)
          message(paste("corresponding to k = ", k," for standard gaussian distribution, recommanded is k in [1,2] for robust results",sep=""))
          message(paste("transferred to tau = ",tau,sep=""))
          
        } else {stop("pi is an amount, choose value in [0,1]")}
      } else {tau=0.8213748
      message(paste("using tau = ",tau,sep=""))
      message("corresponding to k = 1.345")
      message("corresponding to pi = 0.05791489")}
    }
  }
  
  
  Family(
    loss=function(y,f, w = rep(1, length(y))){
      d_1 <- quantile(abs(y-f)[rep(1:length(y), w)], probs = tau)
      ifelse( y-f < -d_1, -2*(1 - q)*d_1*(y - f)-(1 - q)*d_1^2, 
              ifelse( y-f <0,(1 - q)*(y - f)^2, 
                      ifelse(y - f <= d_1, q*(y - f)^2, 2*q*d_1*(y - f)-q*d_1^2)))},
    ngradient = function(y, f, w = rep(1, length(y))){
      d_1 <- quantile(abs(y-f)[rep(1:length(y), w)], probs = tau)
      ifelse( y - f< -d_1, -2*(1 - q)*d_1, ifelse(y - f < 0, 2*(1 - q)*(y - f), ifelse(y - f <= d_1,2*q*(y - f),2*q*d_1)))},
    offset = function(y,w=rep(1, length(y))){
      median(y[rep(1:length(y), w)])}
  )
}

###########################################################################################

# Bisquare
AdaptBisquare <- function (tau = 0.99 , k=NULL, pi=NULL){ # 95% efficiency: k=4.685 translated to tau=0.9999972
  k_to_tau <- function(k){
    return((pnorm(k)-pnorm(-k))) # =1-2*pnorm(-k)
  }
  
  if( !is.null(tau)){
    if((0<=tau)*(tau<=1)){message(paste("using tau = ",tau,sep=""))
      tau_to_k <- function(tau){
        to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
        return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
      }
      k=tau_to_k(tau)
      message(paste("corresponding to k = ",k," for standard gaussian distribution",sep=""))
      pi=1-tau
      message(paste("corresponding to pi = ", pi,sep=""))
    } else {stop("tau is a quantile, choose value in [0,1]")}
  } else{
    if(!is.null(k)){
      if(k>0){message(paste("using k = ",k," for standard gaussian distribution",sep=""))
        tau=k_to_tau(k)
        pi=1-tau
        message(paste("corresponding to pi = ",pi,sep=""))
        message(paste("transferred to tau = ",tau, sep=""))
      } else {stop("k has to be a positive value")}
    } else{
      if(!is.null(pi)){
        if((0<=pi)*(pi<=1)){  message(paste("using pi = ",pi,sep=""))
          tau=1-pi
          tau_to_k <- function(tau){
            to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
            return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
          }
          k=tau_to_k(tau)
          message(paste("corresponding to k = ", k," for standard gaussian distribution",sep=""))
          message(paste("transferred to tau = ",tau,sep=""))
        } else {stop("pi is an amount, choose value in [0,1]")}
      } 
    }
  }
  
  Family(ngradient = function(y, f,w = rep(1, length(y))) {
    d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
    (abs(y-f)<=d)*(y-f)*(1-((y-f)/d)^2)^2}, 
    loss = function(y, f,w = rep(1, length(y))) {
      d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
      (abs(y-f)<=d)*d^2/6*(1-(1-((y-f)/d)^2)^3)+(abs(y-f)>d)/6*1},
    offset = function(y,w=rep(1, length(y))){
      median(y[rep(1:length(y), w)])}
  )
}

###########################################################################################

# Welsh:
AdaptWelsh <- function (tau = 0.99 , k=NULL, pi=NULL){ #95% efficiency: k=2.9846 translated to tau= 0.9971605
  k_to_tau <- function(k){
    return((pnorm(k)-pnorm(-k))) # =1-2*pnorm(-k)
  }
  
  if( !is.null(tau)){
    if((0<=tau)*(tau<=1)){message(paste("using tau = ",tau,sep=""))
      tau_to_k <- function(tau){
        to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
        return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
      }
      k=tau_to_k(tau)
      message(paste("corresponding to k = ",k," for standard gaussian distribution",sep=""))
      pi=1-tau
      message(paste("corresponding to pi = ", pi,sep=""))
    } else {stop("tau is a quantile, choose value in [0,1]")}
  } else{
    if(!is.null(k)){
      if(k>0){message(paste("using k = ",k," for standard gaussian distribution",sep=""))
        tau=k_to_tau(k)
        pi=1-tau
        message(paste("corresponding to pi = ",pi,sep=""))
        message(paste("transferred to tau = ",tau, sep=""))
      } else {stop("k has to be a positive value")}
    } else{
      if(!is.null(pi)){
        if((0<=pi)*(pi<=1)){  message(paste("using pi = ",pi,sep=""))
          tau=1-pi
          tau_to_k <- function(tau){
            to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
            return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
          }
          k=tau_to_k(tau)
          message(paste("corresponding to k = ", k," for standard gaussian distribution",sep=""))
          message(paste("transferred to tau = ",tau,sep=""))
        } else {stop("pi is an amount, choose value in [0,1]")}
      } 
    }
  }
  
  Family(ngradient = function(y, f,w = rep(1, length(y))) {
    d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
    (y-f)*exp(-(y-f)^2/d^2/2)}, 
    loss = function(y, f,w = rep(1, length(y))) {
      d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
      1-exp(-(y-f)^2/d^2/2) },# here we reskaled it if you want to avoid this devide the loss by d^2
    offset = function(y,w=rep(1, length(y))){
      median(y[rep(1:length(y), w)])}
  )
}

###########################################################################################

# Nonnegative Garrote:
AdaptNonnegativeGarrote <- function (tau =0.9544997 , k=NULL, pi=NULL){ #95% efficiency: k=2 translated to tau=0.9544997 
  k_to_tau <- function(k){
    return((pnorm(k)-pnorm(-k))) # =1-2*pnorm(-k)
  }
  
  if( !is.null(tau)){
    if((0<tau)*(tau<=1)){message(paste("using tau = ",tau,sep=""))
      tau_to_k <- function(tau){
        to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
        return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
      }
      k=tau_to_k(tau)
      message(paste("corresponding to k = ",k," for standard gaussian distribution",sep=""))
      pi=1-tau
      message(paste("corresponding to pi = ", pi,sep=""))
    } else {stop("tau is a quantile, choose value in [0,1]")}
  } else{
    if(!is.null(k)){
      if(k>0){message(paste("using k = ",k," for standard gaussian distribution",sep=""))
        tau=k_to_tau(k)
        pi=1-tau
        message(paste("corresponding to pi = ",pi,sep=""))
        message(paste("transferred to tau = ",tau, sep=""))
      } else {stop("k has to be a positive value")}
    } else{
      if(!is.null(pi)){
        if((0<=pi)*(pi<=1)){  message(paste("using pi = ",pi,sep=""))
          tau=1-pi
          tau_to_k <- function(tau){
            to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
            return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
          }
          k=tau_to_k(tau)
          message(paste("corresponding to k = ", k," for standard gaussian distribution",sep=""))
          message(paste("transferred to tau = ",tau,sep=""))
        } else {stop("pi is an amount, choose value in [0,1]")}
      } 
    }
  }
  
  Family(ngradient = function(y, f,w = rep(1, length(y))) {
    d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
    ifelse(abs(y-f)<=d,(y-f),d^2/sign(y-f)/abs(y-f))},
    loss = function(y, f,w = rep(1, length(y))) {
      d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
      ifelse(abs(y-f)<=d,(y-f)^2/2,0.5*d^2*(1+2*log(abs(y-f)/d)))},
    offset = function(y,w=rep(1, length(y))){
      median(y[rep(1:length(y), w)])}
  )
}

###########################################################################################

# Cauchy:
AdaptCauchy <- function (tau = 0.9829162 , k=NULL, pi=NULL){ #95% efficiency: k=2.3849 translated to tau=0.9829162
  k_to_tau <- function(k){
    return((pnorm(k)-pnorm(-k))) # =1-2*pnorm(-k)
  }
  
  if( !is.null(tau)){
    if((0<tau)*(tau<=1)){message(paste("using tau = ",tau,sep=""))
      tau_to_k <- function(tau){
        to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
        return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
      }
      k=tau_to_k(tau)
      message(paste("corresponding to k = ",k," for standard gaussian distribution",sep=""))
      pi=1-tau
      message(paste("corresponding to pi = ", pi,sep=""))
    } else {stop("tau is a quantile, choose value in [0,1]")}
  } else{
    if(!is.null(k)){
      if(k>0){message(paste("using k = ",k," for standard gaussian distribution",sep=""))
        tau=k_to_tau(k)
        pi=1-tau
        message(paste("corresponding to pi = ",pi,sep=""))
        message(paste("transferred to tau = ",tau, sep=""))
      } else {stop("k has to be a positive value")}
    } else{
      if(!is.null(pi)){
        if((0<=pi)*(pi<=1)){  message(paste("using pi = ",pi,sep=""))
          tau=1-pi
          tau_to_k <- function(tau){
            to_minimize <- function(k){return(abs(pnorm(k)-pnorm(-k)-tau))}
            return(optimize(f=to_minimize,interval = c(0.1,8.0),maximum=FALSE)$minimum)
          }
          k=tau_to_k(tau)
          message(paste("corresponding to k = ", k," for standard gaussian distribution",sep=""))
          message(paste("transferred to tau = ",tau,sep=""))
        } else {stop("pi is an amount, choose value in [0,1]")}
      } 
    }
  }
  
  Family(ngradient = function(y, f,w = rep(1, length(y))) {
    d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
    (y-f)/(1+(y-f)^2/d^2)}, 
    loss = function(y, f,w = rep(1, length(y))) {
      d <- quantile(abs(y - f)[rep(1:length(y), w)],probs=tau)
      d^2/2*log(1+(y-f)^2/d^2)},# 
    offset = function(y,w=rep(1, length(y))){
      median(y[rep(1:length(y), w)])}
  )
}
###########################################################################################
