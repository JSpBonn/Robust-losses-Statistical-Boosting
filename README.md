# Robust-losses-Statistical-Boosting

In this Github repository you can find implementations of different quantile-based adaptive loss functions, which are more robust than the L2 loss, for statistical boosting.

Using the R add-on package "mboost", you can specify different loss functions via "family = NEWFAMILY()" as an additional argument inside of the "glmboost" function. For the family itself, there are different arguments, where the user can choose from. A "tau" explains the quantile of the residuals, which are considered non-corrupted, "pi" is a specification for the corrupted amount and "k" was usually used for M-regression, which is translated into "tau" for high asymptotic efficiency, assuming there are Gaussian distributed errors (this is explained in more detail in the article "Robust statistical boosting with quantile-based adaptive loss functions" from Speller et al.).

Example for linear regression:
model <- glmboost(y ~ x1 + x2, data= data.frame(y,x1,x2), family=AdaptHuber(tau=0.8))

Overview adaptive loss functions:
AdaptHuber, AdaptBisquare, AdaptWelsh, AdaptNonnegativeGarrote, AdaptCauchy

E.g. for non-adaptive loss functions:
Gaussian (default for glmboost), Laplace


```{r}
# R code for running an example on a bodyfat data set with different covariates:

source{Robust_quantile_based_adaptive_loss_functions_for_statistical_boosting.R} # loading adaptive loss functions and R-Package "mboost"
data("bodyfat", package = "TH.data") # load bodyfat data

# outcome: DEXfat
# all available covariates (age + waistcirc + hipcirc + elbowbreadth + kneebreadth + anthro3a + anthro3b + anthro3c + anthro4)
set.seed(321)
glm_Huber <- glmboost(DEXfat~. , family = AdaptHuber(tau=0.80) , data = bodyfat)
cvr_H <- cvrisk(glm_Huber,grid = 1:200) # default method is 25-fold bootstrap
coef(glm_Huber[mstop(cvr_H)] , off2int=TRUE , which="") #  coefficients of glm_Huber at optimal stopping iteration for cvrisk

set.seed(321)
glm_Bisquare <- glmboost(DEXfat~. , family = AdaptBisquare(tau=0.99) , data=bodyfat)
cvr_B <- cvrisk(glm_Bisquare,grid = 1:200) # default method is 25-fold bootstrap
coef(glm_Bisquare[mstop(cvr_B)] , off2int=TRUE , which="")  # coefficients of glm_Bisquare at optimal stopping iteration for cvrisk

plot(glm_Huber,ylim=c(-2,5))
plot(glm_Bisquare,ylim=c(-2,5))


# for specific covariates:
set.seed(321)
glm_Huber_2 <- glmboost(DEXfat ~ age + waistcirc + hipcirc , family = AdaptHuber(tau=0.8) , data = bodyfat) 
#cvr_H_2 <- cvrisk(glm_Huber_2,grid = 1:500) # default method is 25-fold bootstrap
#coef(glm_Huber_2[mstop(cvr_H_2)],off2int=TRUE)

set.seed(321)
glm_Bisquare_2 <- glmboost(DEXfat ~ age + waistcirc + hipcirc , family = AdaptBisquare(tau=0.99) , data=bodyfat)
#cvr_B_2 <- cvrisk(glm_Bisquare_2,grid = 1:500) # default method is 25-fold bootstrap
#coef(glm_Bisquare_2[mstop(cvr_B_2)],off2int=TRUE)

```
