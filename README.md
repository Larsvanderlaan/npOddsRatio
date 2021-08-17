# npOddsRatio
Semiparametric Targeted Maximum-Likelihood estimation of the conditional odds ratio in a partially-linear logistic-link regression model with post-treatment informative missingness.

Utilizing the framework of Targeted Maximum-Likelihood estimation (TMLE), efficient estimates and inference is given for the conditional odds ratio in a partially-linear logistic-link regression model with post-treatment informative missingness.
The user can supply an arbitrary parametric form of the conditional log odds ratio. Estimates and confidence intervals are returned. 
Nuisance functions can be estimated using machine-learning (specifically the machine-learning pipeline R package tlverse/sl3), thereby avoiding biases due to model misspecification.
This package allows for the outcome missingness to informed by pre-treatment variables W and the treatment (See function npOR), or through pre-treatment variables, treatment, and post-treatment variables (see function npORMissing).
To estimate the nuisance functions, for convenience, we allow the user to specify either an sl3 Learner object or a R formula object in which case the nuisance functions are parametrically estimated using glm (not recommended). By default, tlverse/hal9001 is used to estimate all nuisances.

Consider the data-structure $(W,A,Y)$ where $W$ is a vector of pre-treatment baseline covariates, $A$ is a binary treatment assignment, and $Y$ is a binary outcome.
We assume the partially-linear logistic-link model
$$\text{logit}(P(Y=1|A=a,W=w)) = a\beta^T\underline{f}(w) + h(w)$$
where $\underline{f}(w)$ is a known vector valued function (the parametric component) and $h(w) = \text{logit}(P(Y=1|A=0,W=w))$ is unknown and unspecified (the nonparametric component).
Note, we have that the conditional odds ratio satisfies
$$\frac{P(Y=1|A=1,W=w)/P(Y=0|A=1,W=w)}{P(Y=1|A=0,W=w)/P(Y=0|A=0,W=w)} = \exp\left\{\beta^T\underline{f}(w) \right\}$$

We are interested in estimating and obtaining inference for the vector $\beta$ and the conditional log odds ratio $\beta^T\underline{f}(w) $. This involves estimating $\beta$ and $h(w)$ first using machine-learning (specifically the 
Highly Adaptive Lasso as implemented in the R package tlverse/hal9001), and then debiasing the initial estimator with the Targeted Maximum-Likelihood estimation (TMLE) framework.
This is implemented in the function "npOR".


In some cases, the outcome $Y$ may be missing. That is, we actually observe $\Delta Y$ where $\Delta$ is a binary variable that is $1$ if $Y$ is observed. 
The missingness $\Delta$ may be informed by both W and A. The conditional odds ratio can still be estimated in this case, and this method is implemented in the function "npOR". 
In more complex settings, one might have that the missingness is informed by a post-treatment variable $Z$. The conditional odds ratio in the presence of post-treatment informed missingness can be etimated using the function "npORMissing".
