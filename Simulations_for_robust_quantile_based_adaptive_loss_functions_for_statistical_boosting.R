# Simulations on "Robust statistical boosting with quantile-based adaptive loss functions:


library(MASS)
library(mvtnorm)
library(mboost)
library(parallel)

source(Robust_quantile_based_adaptive_loss_functions_for_statistical_boosting.R) # loading the adaptive loss functions

###############################################################################
# translation of tau, k, pi:

# translation: k to pi
kpi_func <- function(k){
  return((1/(1-2*pnorm(-k)+2*dnorm(k)/k)-1)*(-1))
}
# translation: k to tau
ktau_func <- function(k){
  return((pnorm(k)-pnorm(-k))) # =1-2*pnorm(-k)
}





###################################################################################
###################################################################################
###################################################################################
Simfunc <- function(id,p=10,n=200){
  # at least p=6 due to number of informative variables
  
  corrupted = c(0.0,0.01,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5) # amount of corruption
  mu = 0 
  sigma = 1 
  correl = 0.5 # correlation of all covariables
  cor_strong = 4 # how large is the corruption? (4 sigma)
  informative_beta = c(1.5,1,0.5,-0.5,-1,-1.5)
  m_stop =3000 # maximal stopping iteration
  n_2=1000 # number of test observations to validate results
  
  corrupted2 <- c(0.0,0.075,0.15) # current different amounts of corruption
  high_dimensional <- p > n # 
  numbercorrupted <- ceiling(n*corrupted2) # number of corrupted observations is an integer
  numberofmodels <-  1+1+2+ 1+3+ (1+1+1)*(1-high_dimensional) # number of models (some deleted due to testing)  
  
  out.tab <- vector("list",14) # one list for every level of "corrupted", not necessary anymore to have 14 
  
  
  tau2 <- c(0.8213748,ktau_func(1),ktau_func(2))# settings Huber
  tau2bis <-c(0.99,0.95,0.90) # settings bisquare
  
  for (i in 1:length(corrupted2)) { # loop for different amount of corruption in the data
    
    set.seed(i*1000+id+1000)
    # simulation setting, run number is id:
    t=length(informative_beta)
    beta <- c(informative_beta,rep(0,p-t))
    
    mean_X = rep(mu,length = p)
    sigma_X = matrix(correl,p,p)
    diag(sigma_X) = sigma
    
    # trainingsdata
    X = rmvnorm(n,mean = mean_X, sigma = sigma_X)
    Summe  = X %*% beta
    sigma_e =sigma
    y = rnorm(n,Summe, sigma_e)
    sigma_ver <- 4
    if(numbercorrupted[i]>0) y[1:numbercorrupted[i]] <-  y[1:numbercorrupted[i]]+rnorm(numbercorrupted[i],0,sigma_ver) ### YYY non-systematic corruption, comment for systematic corruption, choose only one type of corruption
    # if(numbercorrupted[i]>0) y[1:numbercorrupted[i]] <- y[1:numbercorrupted[i]]+cor_strong*sd(y)*(-1)^rbinom(numbercorrupted[i], size=1, prob=0.5) ### YYY systematic corruption, uncomment for systematic corruption, choose only one type of corruption
    datasimulation <- data.frame(y,X)
    
    
    # testdata
    set.seed(i*1000+50000+id)
    mean_X_test = rep(mu,length = p)
    sigma_X_test = matrix(correl,p,p)
    diag(sigma_X_test) = sigma
    X_test = rmvnorm(n_2,mean = mean_X_test, sigma = sigma_X_test)
    Summe_test  = X_test %*% beta
    sigma_e =sigma
    y_test = rnorm(n_2,Summe_test, sigma_e)
    datasimulation_test <- data.frame(y_test,X_test)
    colnames(datasimulation_test) <- c("y","X")
    names(datasimulation_test) <- names(datasimulation)
    
    
    ####################################################################################################################
    rlm_fit <- rlm(y~.,data=datasimulation,maxit=100) ### XXX comment for high dimensional setting
    lm_fit <- lm(y~.,data=datasimulation) ### XXX comment for high dimensional setting
    rlm_bisquare <- rlm(y~.,method="MM" ,data=datasimulation,maxit=100) ### XXX comment high dimensional setting
    
    coefmatrix_rlm <- matrix(0,ncol=p+1,nrow=1) ### XXX comment for high dimensional setting
    coefmatrix_rlm <- rlm_fit$coefficients ### XXX comment for high dimensional setting
    coefmatrix_lm <- matrix(0,ncol=p+1,nrow=1) ### XXX comment for high dimensional setting
    coefmatrix_lm <- lm_fit$coefficients ### XXX comment high dimensional setting
    coefmatrix_bisquare <- matrix(0,ncol=p+1,nrow=1) ### XXX comment for high dimensional setting
    coefmatrix_bisquare <- rlm_bisquare$coefficients ### XXX comment for high dimensional setting
    ####################################################################################################################
    
    # save also estimates, which are "converged"
    coeffmatrix_converg10khub <- matrix(0,ncol =p+1 ,nrow = 3)
    coeffmatrix_optistophub <- matrix(0,ncol =p+1 ,nrow = 3)
    coeffmatrix_converg10kbis <- matrix(0,ncol =p+1 ,nrow = 3)
    coeffmatrix_optistopbis <- matrix(0,ncol =p+1 ,nrow = 3)
    coeffmatrix_rest <- matrix(0,ncol =p+1 ,nrow = 4)
    
    
    MAErlm <- data.frame(matrix(0,ncol=3,nrow=1))
    MAErlm[1,1] <- mean(abs(predict(rlm_fit,newdata=datasimulation_test)-datasimulation_test$y))
    MAErlm[1,2] <- mean(abs(predict(lm_fit,newdata=datasimulation_test)-datasimulation_test$y))
    MAErlm[1,3] <- mean(abs(predict(rlm_bisquare,newdata=datasimulation_test)-datasimulation_test$y))
    colnames(MAErlm) <- c("MAEhub","MAElm","MAEbisquare")
    
    MAEHub <- data.frame(matrix(0, ncol=3,nrow=1))
    MAEBis <- data.frame(matrix(0, ncol=3,nrow=1))
    colnames(MAEHub) <- c("MAEHub default","MAEHub 1", "MAEHub 2") 
    colnames(MAEBis) <- c("MAEBis default",paste("MAEBis",corrupted2[2:3],sep="_"))
    
    cvr_hub <- 1:3 # index corruption
    cvr_bis <- 1:3 # index corruption
    
    
    
    # models:
    
    
    # Huber default:
    set.seed(i*1000+150000+id)
    tau_huber <- 0.8213748 # default huber 
    boost_Huber             <- glmboost(y~.,family=AdaptHuber(tau=tau_huber),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_Huber         <- cvrisk(boost_Huber,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_Huber)>m_stop-50) cvr_boost_Huber <- cvrisk(boost_Huber,grid = 1:8000,papply="lapply")
    
    MAEHub[1,1] <- mean(abs(predict(boost_Huber[mstop(cvr_boost_Huber)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistophub[1,] <- coef(boost_Huber[mstop(cvr_boost_Huber)],which="",off2int=TRUE)
    coeffmatrix_converg10khub[1,] <- coef(boost_Huber[10000],which="",off2int=TRUE) # mstop 10000 saved
    cvr_hub[1] <- mstop(cvr_boost_Huber)
    ###########################
    
    # Huber k=1
    set.seed(i*1000+150000+id)
    tau_huber <-tau2[2]
    boost_Huber             <- glmboost(y~.,family=AdaptHuber(tau=tau_huber),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_Huber         <- cvrisk(boost_Huber,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_Huber)>m_stop-50) cvr_boost_Huber <- cvrisk(boost_Huber,grid = 1:8000,papply="lapply")
    
    MAEHub[1,2] <- mean(abs(predict(boost_Huber[mstop(cvr_boost_Huber)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistophub[2,] <- coef(boost_Huber[mstop(cvr_boost_Huber)],which="",off2int=TRUE)
    coeffmatrix_converg10khub[2,] <- coef(boost_Huber[10000],which="",off2int=TRUE) # mstop 10000 saved
    cvr_hub[2] <- mstop(cvr_boost_Huber)
    ###########################
    
    # Huber k=2
    set.seed(i*1000+150000+id)
    tau_huber <-tau2[3]    
    boost_Huber             <- glmboost(y~.,family=AdaptHuber(tau=tau_huber),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_Huber         <- cvrisk(boost_Huber,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_Huber)>m_stop-50) cvr_boost_Huber <- cvrisk(boost_Huber,grid = 1:8000,papply="lapply")
    
    MAEHub[1,3] <- mean(abs(predict(boost_Huber[mstop(cvr_boost_Huber)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistophub[3,] <- coef(boost_Huber[mstop(cvr_boost_Huber)],which="",off2int=TRUE)
    coeffmatrix_converg10khub[3,] <- coef(boost_Huber[10000],which="",off2int=TRUE) # mstop 10000 saved
    cvr_hub[3] <- mstop(cvr_boost_Huber)
    ###########################
    
    
    
    # bisquare default:
    set.seed(i*1000+250000+id)
    tau_bisquare <- 0.99
    
    boost_bisquare        <- glmboost(y~.,family=AdaptBisquare(tau=tau_bisquare),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_bisquare     <- cvrisk(boost_bisquare,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_bisquare)>m_stop-50) cvr_boost_bisquare <- cvrisk(boost_bisquare,grid = 1:8000,papply="lapply")
    
    MAEBis[1,1] <- mean(abs(predict(boost_bisquare[mstop(cvr_boost_bisquare)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistopbis[1,] <- coef(boost_bisquare[mstop(cvr_boost_bisquare)],which="",off2int=TRUE)
    coeffmatrix_converg10kbis[1,]<- coef(boost_bisquare[10000],which="",off2int=TRUE) # mstop 10000 saved
    
    cvr_bis[1] <- mstop(cvr_boost_bisquare)
    
    ##############################
    set.seed(i*1000+250000+id)
    tau_bisquare <- 0.95
    boost_bisquare        <- glmboost(y~.,family=AdaptBisquare(tau=tau_bisquare),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_bisquare     <- cvrisk(boost_bisquare,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_bisquare)>m_stop-50) cvr_boost_bisquare <- cvrisk(boost_bisquare,grid = 1:8000,papply="lapply")
    
    MAEBis[1,2] <- mean(abs(predict(boost_bisquare[mstop(cvr_boost_bisquare)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistopbis[2,] <- coef(boost_bisquare[mstop(cvr_boost_bisquare)],which="",off2int=TRUE)
    coeffmatrix_converg10kbis[2,]<- coef(boost_bisquare[10000],which="",off2int=TRUE) # mstop 10000 saved
    
    cvr_bis[2] <- mstop(cvr_boost_bisquare)
    
    ###############################
    set.seed(i*1000+250000+id)
    tau_bisquare <- 0.9
    boost_bisquare        <- glmboost(y~.,family=AdaptBisquare(tau=tau_bisquare),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_bisquare     <- cvrisk(boost_bisquare,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_bisquare)>m_stop-50) cvr_boost_bisquare <- cvrisk(boost_bisquare,grid = 1:8000,papply="lapply")
    
    MAEBis[1,3] <- mean(abs(predict(boost_bisquare[mstop(cvr_boost_bisquare)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistopbis[3,] <- coef(boost_bisquare[mstop(cvr_boost_bisquare)],which="",off2int=TRUE)
    coeffmatrix_converg10kbis[3,]<- coef(boost_bisquare[10000],which="",off2int=TRUE) # mstop 10000 saved
    
    cvr_bis[3] <- mstop(cvr_boost_bisquare)
    
    ###############################
    # Laplace and Gaussian
    set.seed(i*1000+250000+id+10400)
    boost_laplace           <- glmboost(y~.,family=Laplace(),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_laplace       <- cvrisk(boost_laplace,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_laplace)>m_stop-50) cvr_boost_laplace <- cvrisk(boost_laplace,grid = 1:8000,papply="lapply")
    
    set.seed(i*1000+250000+id+10600)
    boost_ls                <- glmboost(y~.,control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_ls            <- cvrisk(boost_ls,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_ls)>m_stop-50) cvr_boost_ls <- cvrisk(boost_ls,grid = 1:8000,papply="lapply")
    
    # save estimated coefficients
    coeffmatrix_rest[3,] <- coef(boost_laplace[mstop(cvr_boost_laplace)],which="",off2int=TRUE)
    coeffmatrix_rest[4,] <- coef(boost_ls[mstop(cvr_boost_ls)],which="",off2int=TRUE)
    
    # save MAE
    MAERest <- data.frame(matrix(0,nrow=1, ncol=4))
  
    MAERest[1,3] <- mean(abs(predict(boost_laplace[mstop(cvr_boost_laplace)],newdata=datasimulation_test)-datasimulation_test$y))
    MAERest[1,4] <- mean(abs(predict(boost_ls[mstop(cvr_boost_ls)],newdata=datasimulation_test)-datasimulation_test$y))
    colnames(MAERest) <- c("empty","empty2","MAE_lapl","MAE_ls")
    
    cvr_rest <- c(0,0,mstop(cvr_boost_laplace),mstop(cvr_boost_ls))
    
    # save results depending on loop over i, the index of the amount of corruption in the data
    out.tab[[i]] <- list(corrupted = corrupted[i],id=id,
                         coefmatrix_rlm=coefmatrix_rlm, ### XXX comment for high dimensional setting
                         coefmatrix_lm=coefmatrix_lm, ### XXX comment for high dimensional setting
                         coefmatrix_bisquare=coefmatrix_bisquare, ### XXX comment for high dimensional setting
                         coeffmatrix_converg10khub=coeffmatrix_converg10khub,
                         coeffmatrix_optistophub=coeffmatrix_optistophub,
                         coeffmatrix_converg10kbis=coeffmatrix_converg10kbis,
                         coeffmatrix_optistopbis=coeffmatrix_optistopbis,
                         coeffmatrix_rest=coeffmatrix_rest,
                         cvr_hub=cvr_hub,
                         cvr_bis=cvr_bis,
                         cvr_rest=cvr_rest,
                         MAErlm=MAErlm, ### XXX comment for high dimensional setting
                         MAEHub=MAEHub,
                         MAEBis=MAEBis,
                         MAERest=MAERest
                         ) 
    
  }# end of loop over i
  save(out.tab=out.tab,file = file.path("//home/XXX.../50to200",paste(id,"simulationsERG", "par",p,"observ",n, "sim.RData",sep="_")))
  # save(out.tab=out.tab,file = file.path("//home/XXX.../500to400",paste(id,"simulationsERG", "par",p,"observ",n, "sim.RData",sep="_"))) ### XXX uncomment for high dimensional setting
  return(out.tab)
  # every out.tab is for one id of the loop
}
##########################################################################################################
##########################################################################################################

# chose simulation setting:
# low dimensional:
p =50 ### XXX comment for high dimensional setting
n=200 ### XXX comment for high dimensional setting

# high dimensional:
# take care: some matrices are not going to work for this case (classical regression estimates), look for "XXX uncomment/comment for high dimensional setting" in the code above
# p=500 ### XXX uncomment for high dimensional setting
# n=400 ### XXX uncomment for high dimensional setting

id <-1:100 # 100 runs, id is specific run 



# choose low or high dimensional setting (uncomment/comment XXX) and choose systematic or unsystematic corruption (uncomment/comment YYY) 
ERG <- mclapply(id,FUN=Simfunc,n=n,p=p,mc.cores =25,mc.set.seed = TRUE,  mc.preschedule = FALSE) # parallel setting for cluster (here 25 cores used)


save(ERG,file=paste("simulationsERG", "par",p,"runs",max(id),"observ",n, "sim.RData",sep="_"))



