# Simulations on "Robust statistical boosting with quantile-based adaptive loss functions:


library(MASS)
library(mvtnorm)
library(mboost)
library(parallel)

library(quantreg)
library(glmnet)

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
  
  corrupted = c(0.0,0.075,0.15) # amount of corruption
  mu = 0 
  sigma = 1 
  correl = 0.5 # correlation of all covariables, uncomment for Toeplitz
  # Toeplitz:
  # corr=0.8 # Toeplitz correlation between neighboured predictor variables, comment for Toeplitz
  # help = numeric(p) comment for Toeplitz
  # for (k in 1:p){ help[k]=corr^(k-1) comment for Toeplitz
  # } comment for Toeplitz
  # toeplitzcor = stats::toeplitz(help) comment for Toeplitz
  
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
    sigma_X = matrix(correl,p,p) # comment for Toeplitz
    # sigma_X = matrix(toeplitzcor,p,p) # uncomment for Toeplitz
    diag(sigma_X) = sigma
    
    # trainingsdata
    X = rmvnorm(n,mean = mean_X, sigma = sigma_X)
    Summe  = X %*% beta
    sigma_e =sigma
    y = rnorm(n,Summe, sigma_e)
    sigma_ver <- 4  
       
    
    # symmetric corruption:
    if(numbercorrupted[i]>0) y[1:numbercorrupted[i]] <-  y[1:numbercorrupted[i]]+rnorm(numbercorrupted[i],0,sigma_ver) ### YYY non-systematic corruption, comment for systematic corruption, choose only one type of corruption
    # if(numbercorrupted[i]>0) y[1:numbercorrupted[i]] <- y[1:numbercorrupted[i]]+cor_strong*sd(y)*(-1)^rbinom(numbercorrupted[i], size=1, prob=0.5) ### YYY systematic corruption, uncomment for systematic corruption, choose only one type of corruption
    
    # non-symmetric corruption:
    # if(numbercorrupted[i]>0) y[1:numbercorrupted[i]] <-  y[1:numbercorrupted[i]]+ abs(rnorm(numbercorrupted[i],0,sigma_ver)) ### YYY non-symmetric, non-systematic corruption, comment for systematic corruption, choose only one type of corruption
    # if(numbercorrupted[i]>0) y[1:numbercorrupted[i]] <- y[1:numbercorrupted[i]]+cor_strong*sd(y) ### YYY non-symmetric, systematic corruption, comment for non-systematic corruption, choose only one type of corruption
    
    datasimulation <- data.frame(y,X)
    
    
    # testdata
    set.seed(i*1000+50000+id)
    mean_X_test = rep(mu,length = p)
    sigma_X_test = matrix(correl,p,p) # comment for Toeplitz
    # sigma_X_test = toeplitzcor # uncomment for Toeplitz
    diag(sigma_X_test) = sigma
    X_test = rmvnorm(n_2,mean = mean_X_test, sigma = sigma_X_test)
    Summe_test  = X_test %*% beta
    sigma_e =sigma
    y_test = rnorm(n_2,Summe_test, sigma_e)
    datasimulation_test <- data.frame(y_test,X_test)
    colnames(datasimulation_test) <- c("y","X")
    names(datasimulation_test) <- names(datasimulation)
    
    
    
    # number_methods <-13 # running time
    # time_matrix <- matrix(0,ncol = 2,nrow = number_methods) 
    # time_matrix[,1] <- c("rlm_huber","lm","rlm_bisquare","Huber_def","Huber_k1","Huberk2","Bisquare0_99","Bisquare0_95","Bisquare0_90","L1_boost","L2_boost","lassoquantreg","lasso")
    
    
    ####################################################################################################################
    # time_a <-proc.time()[3] 
    rlm_fit <- rlm(y~.,data=datasimulation,maxit=100) ### XXX comment for (ultra) high dimensional setting
    # time_b <- proc.time()[3]
    # time_matrix[1,2] <- time_b-time_a  
    
    # time_a <-proc.time()[3]
    lm_fit <- lm(y~.,data=datasimulation) ### XXX comment for (ultra) high dimensional setting
    # time_b <- proc.time()[3] 
    # time_matrix[2,2] <- time_b-time_a
    
    # time_a <-proc.time()[3]
    rlm_bisquare <- rlm(y~.,method="MM" ,data=datasimulation,maxit=100) ### XXX (ultra) comment high dimensional setting
    # time_matrix[3,2] <- time_b-time_a
    
    coefmatrix_rlm <- matrix(0,ncol=p+1,nrow=1) ### XXX comment for (ultra) high dimensional setting
    coefmatrix_rlm <- rlm_fit$coefficients ### XXX comment for (ultra) high dimensional setting
    coefmatrix_lm <- matrix(0,ncol=p+1,nrow=1) ### XXX comment for (ultra) high dimensional setting
    coefmatrix_lm <- lm_fit$coefficients ### XXX comment high dimensional setting
    coefmatrix_bisquare <- matrix(0,ncol=p+1,nrow=1) ### XXX comment for (ultra) high dimensional setting
    coefmatrix_bisquare <- rlm_bisquare$coefficients ### XXX comment for (ultra) high dimensional setting
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
    # time_a <-proc.time()[3]
    boost_Huber             <- glmboost(y~.,family=AdaptHuber(tau=tau_huber),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_Huber         <- cvrisk(boost_Huber,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_Huber)>m_stop-50) cvr_boost_Huber <- cvrisk(boost_Huber,grid = 1:8000,papply="lapply")
    # time_b <- proc.time()[3]
    # time_matrix[4,2] <- time_b-time_a 
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
    # time_b <- proc.time()[3]
    # time_matrix[5,2] <- time_b-time_a
    MAEHub[1,2] <- mean(abs(predict(boost_Huber[mstop(cvr_boost_Huber)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistophub[2,] <- coef(boost_Huber[mstop(cvr_boost_Huber)],which="",off2int=TRUE)
    coeffmatrix_converg10khub[2,] <- coef(boost_Huber[10000],which="",off2int=TRUE) # mstop 10000 saved
    cvr_hub[2] <- mstop(cvr_boost_Huber)
    ###########################
    
    # Huber k=2
    set.seed(i*1000+150000+id)
    tau_huber <-tau2[3]
    # time_a <-proc.time()[3]
    boost_Huber             <- glmboost(y~.,family=AdaptHuber(tau=tau_huber),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_Huber         <- cvrisk(boost_Huber,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_Huber)>m_stop-50) cvr_boost_Huber <- cvrisk(boost_Huber,grid = 1:8000,papply="lapply")
    # time_b <- proc.time()[3]
    # time_matrix[6,2] <- time_b-time_a
    MAEHub[1,3] <- mean(abs(predict(boost_Huber[mstop(cvr_boost_Huber)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistophub[3,] <- coef(boost_Huber[mstop(cvr_boost_Huber)],which="",off2int=TRUE)
    coeffmatrix_converg10khub[3,] <- coef(boost_Huber[10000],which="",off2int=TRUE) # mstop 10000 saved
    cvr_hub[3] <- mstop(cvr_boost_Huber)
    ###########################
    
    
    
    # bisquare default:
    set.seed(i*1000+250000+id)
    tau_bisquare <- 0.99
    # time_a <-proc.time()[3]
    boost_bisquare        <- glmboost(y~.,family=AdaptBisquare(tau=tau_bisquare),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_bisquare     <- cvrisk(boost_bisquare,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_bisquare)>m_stop-50) cvr_boost_bisquare <- cvrisk(boost_bisquare,grid = 1:8000,papply="lapply")
    # time_b <- proc.time()[3]
    # time_matrix[7,2] <- time_b-time_a
    MAEBis[1,1] <- mean(abs(predict(boost_bisquare[mstop(cvr_boost_bisquare)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistopbis[1,] <- coef(boost_bisquare[mstop(cvr_boost_bisquare)],which="",off2int=TRUE)
    coeffmatrix_converg10kbis[1,]<- coef(boost_bisquare[10000],which="",off2int=TRUE) # mstop 10000 saved
    
    cvr_bis[1] <- mstop(cvr_boost_bisquare)
    
    ##############################
    set.seed(i*1000+250000+id)
    tau_bisquare <- 0.95
    # time_a <-proc.time()[3]
    boost_bisquare        <- glmboost(y~.,family=AdaptBisquare(tau=tau_bisquare),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_bisquare     <- cvrisk(boost_bisquare,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_bisquare)>m_stop-50) cvr_boost_bisquare <- cvrisk(boost_bisquare,grid = 1:8000,papply="lapply")
    # time_b <- proc.time()[3]
    # time_matrix[8,2] <- time_b-time_a
    MAEBis[1,2] <- mean(abs(predict(boost_bisquare[mstop(cvr_boost_bisquare)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistopbis[2,] <- coef(boost_bisquare[mstop(cvr_boost_bisquare)],which="",off2int=TRUE)
    coeffmatrix_converg10kbis[2,]<- coef(boost_bisquare[10000],which="",off2int=TRUE) # mstop 10000 saved
    
    cvr_bis[2] <- mstop(cvr_boost_bisquare)
    
    ###############################
    set.seed(i*1000+250000+id)
    tau_bisquare <- 0.9
    # time_a <-proc.time()[3]
    boost_bisquare        <- glmboost(y~.,family=AdaptBisquare(tau=tau_bisquare),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_bisquare     <- cvrisk(boost_bisquare,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_bisquare)>m_stop-50) cvr_boost_bisquare <- cvrisk(boost_bisquare,grid = 1:8000,papply="lapply")
    # time_b <- proc.time()[3]
    # time_matrix[9,2] <- time_b-time_a
    MAEBis[1,3] <- mean(abs(predict(boost_bisquare[mstop(cvr_boost_bisquare)],newdata=datasimulation_test)-datasimulation_test$y))
    
    coeffmatrix_optistopbis[3,] <- coef(boost_bisquare[mstop(cvr_boost_bisquare)],which="",off2int=TRUE)
    coeffmatrix_converg10kbis[3,]<- coef(boost_bisquare[10000],which="",off2int=TRUE) # mstop 10000 saved
    
    cvr_bis[3] <- mstop(cvr_boost_bisquare)
    
    ###############################
    # Laplace and Gaussian
    set.seed(i*1000+250000+id+10400)
    # time_a <-proc.time()[3]
    boost_laplace           <- glmboost(y~.,family=Laplace(),control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_laplace       <- cvrisk(boost_laplace,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_laplace)>m_stop-50) cvr_boost_laplace <- cvrisk(boost_laplace,grid = 1:8000,papply="lapply")
    # time_b <- proc.time()[3]
    # time_matrix[10,2] <- time_b-time_a
    set.seed(i*1000+250000+id+10600)
    # time_a <-proc.time()[3]
    boost_ls                <- glmboost(y~.,control = boost_control(mstop=1),data=datasimulation)
    cvr_boost_ls            <- cvrisk(boost_ls,grid = 1:m_stop,papply="lapply")
    if(mstop(cvr_boost_ls)>m_stop-50) cvr_boost_ls <- cvrisk(boost_ls,grid = 1:8000,papply="lapply")
    # time_b <- proc.time()[3]
    # time_matrix[11,2] <- time_b-time_a
    # save estimated coefficients
    coeffmatrix_rest[3,] <- coef(boost_laplace[mstop(cvr_boost_laplace)],which="",off2int=TRUE)
    coeffmatrix_rest[4,] <- coef(boost_ls[mstop(cvr_boost_ls)],which="",off2int=TRUE)
    
    # save MAE
    MAERest <- data.frame(matrix(0,nrow=1, ncol=4))
  
    MAERest[1,3] <- mean(abs(predict(boost_laplace[mstop(cvr_boost_laplace)],newdata=datasimulation_test)-datasimulation_test$y))
    MAERest[1,4] <- mean(abs(predict(boost_ls[mstop(cvr_boost_ls)],newdata=datasimulation_test)-datasimulation_test$y))
    colnames(MAERest) <- c("empty","empty2","MAE_lapl","MAE_ls")
    
    cvr_rest <- c(0,0,mstop(cvr_boost_laplace),mstop(cvr_boost_ls))
    
    ################################
    
    set.seed(i*1000+250000+id+108000)
    # lasso robust, not sparse: take care about memory issues for ultra high dimensional
    # time_a <- proc.time()[3]
    f_lasso <- rq(y ~ .,data=datasimulation, method="lasso")
    # time_b <- proc.time()[3]
    # time_matrix[12,2] <- time_b-time_a
    
    set.seed(i*1000+250000+id+208000)
    # lasso (classical)
    # time_a <- proc.time()[3]
    mod_cv <- cv.glmnet(x=as.matrix(datasimulation[,c(2:c(p+1))]), y=datasimulation$y, family="gaussian")
    Lasso_beta = coef(mod_cv, mod_cv$lambda.min) 
    # time_b <- proc.time()[3] 
    # time_matrix[13,2] <- time_b-time_a 
       
    
    lasso_coef_vector <- as.vector(Lasso_beta)
    
    MAE_lasso <- matrix(0,ncol = 1,nrow=2)
    MAE_lasso[1,1] <- mean(abs(predict(f_lasso,newdata=datasimulation_test[,2:c(p+1)])-datasimulation_test$y))
    
    sim_dat <- as.data.frame(matrix(0,ncol=p+1,nrow = n_2))
    sim_dat[,1] <- 1
    sim_dat[,c(2:c(p+1))] <- datasimulation_test[,c(2:c(p+1))]
    
    lasso_vec_error <- rep(0,n_2)
    for (j in 1:n_2) {lasso_vec_error[j] <- sum(sim_dat[j,]*as.vector(Lasso_beta))}
    MAE_lasso[2,1] <- mean(abs(lasso_vec_error -datasimulation_test[,1]))
    
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
                         # time_matrix=time_matrix,
                         lasso_coef_vector=lasso_coef_vector,
                         MAE_lasso=MAE_lasso
                         ) 
    
  } # end of loop over i
  save(out.tab=out.tab,file = file.path("//home/XXX.../50to200",paste(id,"simulationsERG", "par",p,"observ",n, "sim.RData",sep="_")))
  # save(out.tab=out.tab,file = file.path("//home/XXX.../500to400",paste(id,"simulationsERG", "par",p,"observ",n, "sim.RData",sep="_"))) ### XXX uncomment for (ultra) high dimensional setting
  # save(out.tab=out.tab,file = file.path("//home/XXX.../15000to200",paste(id,"simulationsERG", "par",p,"observ",n, "sim.RData",sep="_"))) ### XXX uncomment for (ultra) high dimensional setting
 
  return(out.tab)
  # every out.tab is for one id of the loop
}
##########################################################################################################
##########################################################################################################

# chose simulation setting:
# low dimensional:
p =50 ### XXX comment for (ultra) high dimensional setting
n=200 ### XXX comment for (ultra) high dimensional setting

# high dimensional:
# take care: some matrices are not going to work for this case (classical regression estimates), look for "XXX uncomment/comment for high dimensional setting" in the code above
# p=500 ### XXX uncomment for high dimensional setting
# n=400 ### XXX uncomment for high dimensional setting


# ultra high dimensional:
# take care: some matrices are not going to work for this case (classical regression estimates), look for "XXX uncomment/comment for (ultra) high dimensional setting" in the code above
# p=15000 XXX uncomment for ultra high dimensional setting
# n=200 XXX uncomment for ultra high dimensional setting
# # take care: computational issues in the data generating process, we recommand to generate blocks of 1000x1000 for this setting and also to remove all you do not need to save memory:

# # Adjust/ Add this to your code for ultra high dimensional:
#  p_2 <- p/15
# mean_X = rep(mu,length = p_2)
# sigma_X = matrix(correl,p_2,p_2) # or for Toeplitz: matrix(toeplitzcor,p_2,p_2)
# diag(sigma_X)=sigma
# # trainingsdata
# X1 = rmvnorm(n,mean = mean_X, sigma = sigma_X)
# ...
# X15 = rmvnorm(n,mean = mean_X, sigma = sigma_X)
# X=cbind(X1,...,X15)
# testdata
# mean_X_test = rep(mu,length = p_2)
# sigma_X_test = matrix(correl,p_2,p_2) # or for Toeplitz:  matrix(toeplitzcor,p_2,p_2)
# diag(sigma_X_test) = sigma
# X_test1 = rmvnorm(n_2,mean = mean_X_test, sigma = sigma_X_test)
# ...
# X_test15 = rmvnorm(n_2,mean = mean_X_test, sigma = sigma_X_test)
# X_test=cbind(X_test1,...,X_test15)



id <-1:1000 # 1000 runs, id is specific run 



# choose low or high or ultra high dimensional setting (uncomment/comment XXX) and choose (symmetric/non-symmetric) systematic or unsystematic corruption (uncomment/comment YYY) and choose symmetric correlation or Toeplitz correlation structure (uncomment/comment for Toeplitz)
ERG <- mclapply(id,FUN=Simfunc,n=n,p=p,mc.cores =25,mc.set.seed = TRUE,  mc.preschedule = FALSE) # parallel setting for cluster (here 25 cores used)


save(ERG,file=paste("simulationsERG", "par",p,"runs",max(id),"observ",n, "sim.RData",sep="_"))



