# Robust-losses-Statistical-Boosting

In this Github repository you can find implementations of different quantile-based adaptive loss functions, which are more robust than the L2 loss, for statistical boosting.

Using the R add-on package "mboost", you can specify different loss functions via "family = NEWFAMILY()" as an additional argument inside of the "glmboost" function. For the family itself, there are different arguments, where the user can choose from. A "tau" explains the quantile of the residuals, which are considered non-corrupted, "pi" is a specification for the corrupted amount and "k" was usually used for M-regression, which is translated into "tau" for high asymptotic efficiency, assuming there are Gaussian distributed errors (this is explained in more detail in the article "Robust statistical boosting with quantile-based adaptive loss functions" from Speller et al.).

Example for linear regression:
model <- glmboost(y ~ x1 + x2, data= data.frame(y,x1,x2), family=AdaptHuber(tau=0.8))

Overview adaptive loss functions:
AdaptHuber
AdaptBisquare
AdaptWelsh
AdaptNonnegativeGarrote
AdaptCauchy

E.g. for non-adaptive loss functions:
Gaussian   (default for glmboost)
Laplace
