# Github Application:
#library(parallel)
library(glmnet) # normal lasso using glmnet
library(quantreg) # robust lasso using R function: "rq(, method="lasso")" (corresponding to function "rq.fit.lasso")

#setwd("//folderapplication")
source{Robust_quantile_based_adaptive_loss_functions_for_statistical_boosting.R} # loading adaptive loss functions and R-Package "mboost"

#setwd("//folderapplication")
load(data_set_application.RData) # the loaded data is called: "data_set" 


#########################################################################################
#########################################################################################





number_cancer_cell_lines <- 59
p=14951


#########################################################################################

mstopmatrix <- matrix(0,ncol=8,nrow = number_cancer_cell_lines)
mstopmatrix <- data.frame(mstopmatrix)
colnames(mstopmatrix) <- c("Gaussian","Hub_def","Hub_1","Hub_2","Bis_def","Bis_0.95","Bis_0.90","Laplace")
errormatrix <- matrix(0,ncol=8,nrow = number_cancer_cell_lines)
errormatrix <- data.frame(errormatrix)
colnames(errormatrix) <- c("Gaussian","Hub_def","Hub_1","Hub_2","Bis_def","Bis_0.95","Bis_0.90","Laplace")
time_matrix <- matrix(0,ncol=8,nrow=number_cancer_cell_lines)
colnames(time_matrix) <- c("Gaussian","Hub_def","Hub_1","Hub_2","Bis_def","Bis_0.95","Bis_0.90","Laplace")

# paste("matrix_",c("Gaussian","Hub_def","Hub_1","Hub_2","Bis_def","Bis_0.95","Bis_0.90","Laplace")[1],sep="")
matrix_all<- data.frame(matrix(0,ncol=p+1,nrow=number_cancer_cell_lines))

colnames(matrix_all) <- c("Intercept",colnames(data_set[,-1]))
matrix_Gaussian <- matrix_all
matrix_Hub_def  <- matrix_all
matrix_Hub_1    <- matrix_all
matrix_Hub_2    <- matrix_all
matrix_Bis_def  <- matrix_all
matrix_Bis_0.95 <- matrix_all
matrix_Bis_0.90 <- matrix_all
matrix_Laplace  <- matrix_all


for (i in 1:number_cancer_cell_lines) {
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm1 <- glmboost(KRT19_protein ~ .,data =data_set[-i,],control =boost_control(mstop=1))
  cvr_glm1        <- cvrisk(glm1,grid = 1:200)
  time_b <- proc.time()[3]
  time_matrix[i,1] <- time_b-time_a
  mstopmatrix[i,1] <- mstop(cvr_glm1)
  errormatrix[i,1] <- predict(glm1[mstop(cvr_glm1)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Gaussian[i,] <- coef(glm1[mstop(cvr_glm1)],off2int=TRUE,which="")
  
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm2 <- glmboost( KRT19_protein~.,family=AdaptHuber(),data =data_set[-i,],control=boost_control(mstop=1)) # AdaptHuber: k = 1.345
  cvr_glm2        <- cvrisk(glm2,grid = 1:200)
  time_b <- proc.time()[3]
  time_matrix[i,2] <- time_b-time_a
  mstopmatrix[i,2] <- mstop(cvr_glm2)
  errormatrix[i,2] <- predict(glm2[mstop(cvr_glm2)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Hub_def[i,] <- coef(glm2[mstop(cvr_glm2)],off2int=TRUE,which="")
  
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm3 <- glmboost( KRT19_protein~.,family=AdaptHuber(k=1),data =data_set[-i,],control=boost_control(mstop=1)) # AdaptHuber: k = 1
  cvr_glm3        <- cvrisk(glm3,grid = 1:600)
  time_b <- proc.time()[3]
  time_matrix[i,3] <- time_b-time_a
  mstopmatrix[i,3] <- mstop(cvr_glm3)
  errormatrix[i,3] <- predict(glm3[mstop(cvr_glm3)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Hub_1[i,] <- coef(glm3[mstop(cvr_glm3)],off2int=TRUE,which="")
  
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm4 <- glmboost( KRT19_protein~.,family=AdaptHuber(k=2),data =data_set[-i,],control=boost_control(mstop=1)) # AdaptHuber: k = 2
  cvr_glm4        <- cvrisk(glm4,grid = 1:200)
  time_b <- proc.time()[3]
  time_matrix[i,4] <- time_b-time_a
  mstopmatrix[i,4] <- mstop(cvr_glm4)
  errormatrix[i,4] <- predict(glm4[mstop(cvr_glm4)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Hub_2[i,] <- coef(glm4[mstop(cvr_glm4)],off2int=TRUE,which="")
  
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm5 <- glmboost( KRT19_protein ~ .,family=AdaptBisquare(),data =data_set[-i,],control=boost_control(mstop=1)) # AdaptBisquare: tau = 0.99
  cvr_glm5        <- cvrisk(glm5,grid = 1:200)
  time_b <- proc.time()[3]
  time_matrix[i,5] <- time_b-time_a
  mstopmatrix[i,5] <- mstop(cvr_glm5)
  errormatrix[i,5] <- predict(glm5[mstop(cvr_glm5)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Bis_def[i,] <- coef(glm5[mstop(cvr_glm5)],off2int=TRUE,which="")
  
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm6 <- glmboost( KRT19_protein ~ .,family=AdaptBisquare(0.95),data =data_set[-i,],control=boost_control(mstop=1))  # AdaptBisquare: tau = 0.95
  cvr_glm6        <- cvrisk(glm6,grid = 1:200)
  time_b <- proc.time()[3]
  time_matrix[i,6] <- time_b-time_a
  mstopmatrix[i,6] <- mstop(cvr_glm6)
  errormatrix[i,6] <- predict(glm6[mstop(cvr_glm6)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Bis_0.95[i,] <- coef(glm6[mstop(cvr_glm6)],off2int=TRUE,which="")
  
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm7 <- glmboost( KRT19_protein ~ .,family=AdaptBisquare(0.9),data =data_set[-i,],control=boost_control(mstop=1))  # AdaptBisquare: tau = 0.90
  cvr_glm7        <- cvrisk(glm7,grid = 1:200)
  time_b <- proc.time()[3]
  time_matrix[i,7] <- time_b-time_a
  mstopmatrix[i,7] <- mstop(cvr_glm7)
  errormatrix[i,7] <- predict(glm7[mstop(cvr_glm7)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Bis_0.90[i,] <- coef(glm7[mstop(cvr_glm7)],off2int=TRUE,which="")
  
  set.seed(100+i)
  time_a<-proc.time()[3]
  glm8 <- glmboost( KRT19_protein ~ .,family=Laplace(),data =data_set[-i,],control=boost_control(mstop=1)) # Laplace
  cvr_glm8        <- cvrisk(glm8,grid = 1:200)
  time_b <- proc.time()[3]
  time_matrix[i,8] <- time_b-time_a
  mstopmatrix[i,8] <- mstop(cvr_glm8)
  errormatrix[i,8] <- predict(glm8[mstop(cvr_glm8)],newdata=data_set[i,])-data_set$KRT19_protein[i]
  matrix_Laplace[i,] <- coef(glm8[mstop(cvr_glm8)],off2int=TRUE,which="")
  
  print(i)
}

#setwd("//folderapplication")
#save(time_matrix,file="time_matrix.RData")
#save(errormatrix,file="errormatrix.RData")
#save(mstopmatrix,file="mstopmatrix.RData")

#save(matrix_Gaussian,file="matrix_Gaussian.RData")
#save(matrix_Hub_2,file="matrix_Hub_2.RData")
#save(matrix_Hub_1,file="matrix_Hub_1.RData")
#save(matrix_Hub_def,file="matrix_Hub_def.RData")
#save(matrix_Bis_def,file="matrix_Bis_def.RData")
#save(matrix_Bis_0.95,file="matrix_Bis_0.95.RData")
#save(matrix_Bis_0.90,file="matrix_Bis_0.90.RData")
#save(matrix_Laplace,file="matrix_Laplace.RData")

# # summary of results:
# MAE_methods <- matrix(0,ncol=8)
# MeanMstop <- matrix(0,ncol=8)
# colnames(MAE_methods  ) <- c("Gaussian","Hub_def","Hub_1","Hub_2","Bis_def","Bis_0.95","Bis_0.90","Laplace")
# colnames(MeanMstop) <- c("Gaussian","Hub_def","Hub_1","Hub_2","Bis_def","Bis_0.95","Bis_0.90","Laplace")
# for (j in 1:8) {
#  MAE_methods[1,j] <-  mean(abs(errormatrix[,j]))
#  MeanMstop[1,j] <-  mean(abs(mstopmatrix[,j]))
# }
# 
# 
# View(MAE_methods)
# View(MeanMstop)
# save(MeanMstop,file="Stoppingiteration.RData")
# save(MAE_methods,file="MAE_methods")
# Gaussian_size <- rep(0,59)
# Laplace_size <- rep(0,59)
# Hub_def_size <- rep(0,59) 
# Hub_1_size <-  rep(0,59)
# Hub_2_size <-  rep(0,59)
# Bis_def_size <-  rep(0,59)
# Bis_0.95_size <-  rep(0,59)
# Bis_0.90_size <-  rep(0,59)
# 
# for (i in 1:59) {
#   
#   Gaussian_size[i] <- sum(matrix_Gaussian[i,-1]!=0)
#   Laplace_size[i] <- sum(matrix_Laplace[i,-1]!=0)
#   Hub_def_size[i] <- sum(matrix_Hub_def[i,-1]!=0)
#   Hub_1_size[i] <- sum(matrix_Hub_1[i,-1]!=0)
#   Hub_2_size[i] <- sum(matrix_Hub_2[i,-1]!=0)
#   Bis_def_size[i] <- sum(matrix_Bis_def[i,-1]!=0)
#   Bis_0.95_size[i] <- sum(matrix_Bis_0.95[i,-1]!=0)
#   Bis_0.90_size[i] <- sum(matrix_Bis_0.90[i,-1]!=0)
# }
# 
# c(mean(Gaussian_size), sd(Gaussian_size))
# c(mean(Laplace_size), sd(Laplace_size))
# c(mean(Hub_def_size), sd(Hub_def_size))
# c(mean(Hub_1_size), sd(Hub_1_size))
# c(mean(Hub_2_size), sd(Hub_2_size))
# c(mean(Bis_def_size), sd(Bis_def_size))
# c(mean(Bis_0.95_size), sd(Bis_0.95_size))
# c(mean(Bis_0.90_size), sd(Bis_0.90_size))


#########################################################################################

# Lasso L2:

error_lasso_L2 <- rep(0,number_cancer_cell_lines)
time_lasso_L2 <- rep(0,number_cancer_cell_lines)
number_selected_lasso_L2 <- rep(0,number_cancer_cell_lines)


for (i in 1:number_cancer_cell_lines) {
  set.seed(i+100)
  
  time_a<-proc.time()[3]
  lasso <- cv.glmnet(x=as.matrix(data_set[-i,-1]),y=as.matrix(data_set[-i,1]))
  time_b <- proc.time()[3]
  time_lasso_L2[i] <- time_b-time_a
  
  error_lasso_L2[i]<- predict(lasso,newx = as.matrix(data_set[i,-1]),s="lambda.min")-data_set$KRT19_protein[i]
  number_selected_lasso_L2[i] <- sum(coef(lasso)[2:c(p+1)]!=0)
  print(i)
}
# 
# MAE_lasso_L2 <- mean(abs(error_lasso_L2))
# # MAE_lasso # 32.65
# SD_lasso_L2 <- sd(error_lasso_L2) 
# # SD_lasso_L2 # 44.38285
# mean_time_lasso_L2 <- mean(time_lasso_L2)
# # mean_time_lasso_L2
# SD_time_lasso_L2 <- sd(time_lasso_L2)
# # SD_time_lasso_L2
# mean_number_selected_lasso_L2 <- mean(number_selected_lasso_L2)
# # mean_number_selected_lasso_L2 
# SD_number_selected_lasso_L2 <- sd(number_selected_lasso_L2)
# # SD_number_selected_lasso_L2 


#########################################################################################


# Lasso L1 (robust): 
# take care: there might be memory issues >16GB per robust Lasso model

error_lasso_L1 <- rep(0,number_cancer_cell_lines)
time_lasso_L1 <- rep(0,number_cancer_cell_lines)
coef_matrix <- matrix(0,nrow = number_cancer_cell_lines,ncol = p+1)
number_selected_lasso_L1 <- rep(0,number_cancer_cell_lines)  # all variables get an update



for (i in 1:number_cancer_cell_lines) {
  set.seed(100+i)
  
  time_a<-proc.time()[3]
  f_lasso <- rq(KRT19_protein ~ .,data=as.data.frame(data_set[-i,]), method="lasso")
  time_b <- proc.time()[3]
  time_lasso_L1[i] <- time_b-time_a
  
  coef_matrix[i,] <- f_lasso$coefficients
  error_lasso_L1[i] <- predict(f_lasso,newdata=data_set[i,])-data_set[i,1]
  
  print(i)
}
# 
# MAE_lasso_L1 <- mean(abs(error_lasso_L1))
# # MAE_lasso_L1 # 35.99
# SD_lasso_L1 <- sd(error_lasso_L1)
# # SD_lasso_L1  # 47.81
# mean_time_lasso_L1 <- mean(time_lasso_L1)
# # mean_time_lasso_L1
# SD_time_lasso_L1 <- sd(time_lasso_L1)
# # SD_time_lasso_L1
# 
# # the robust lasso will always estimate all coefficients:
# # number of coefficients < 0.00001:
# for (i in 1:number_cancer_cell_lines) {
#   number_selected_lasso_L1[i] <- sum(abs(coef_matrix[i,2:14952])>0.00001)
# }
# 
# mean_number_selected_lasso_L1 <- mean(number_selected_lasso_L1)
# # mean_number_selected_lasso_L1 
# SD_number_selected_lasso_L1 <- sd(number_selected_lasso_L1)
# # SD_number_selected_lasso_L1 
# 



