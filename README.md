# npOddsRatio
Semiparametric and Nonparametric Targeted Maximum-Likelihood inference for the conditional odds ratio with post-and-pre-treatment informative outcome missingness using black-box machine-learning. 

## Key words: 
Causal Inference, Missing Data, Partially linear logistic regression, semiparametric, nonparametric, TMLE, machine-learning, Targeted-learning, working-model, Highly Adaptive Lasso


## Method implemented
Utilizing the framework of Targeted Maximum-Likelihood estimation (TMLE), efficient estimates and inference is given for the conditional odds ratio in both semiparametric (sp) and nonparametric (np) models with post-treatment informative missingness. Black-box machine-learning is leveraged to estimate nuisance parameters. In short, the methods implemented here can be viewed as (causal) semi/nonparametric generalizations of ordinary parametric logistic regression (e.g. glm). 

## Data structure and target parameter
Consider the data-structure (W,A, Z, Delta, Y) where W is a vector of pre-treatment baseline covariates, A is a binary treatment assignment, and Y is a binary outcome. Delta is a binary variable that captures whether Y is missing (Delta =1 is not missing). Z is a post-treatment confounding variable that informs both Delta (i.e. the missingness of Y) and Y itself. 

The target parameter is the conditional odds ratio:
OR(w): = {P(Y=1|A=1,W=w)/P(Y=0|A=1,W=w)} / {P(Y=1|A=0,W=w)/P(Y=0|A=0,W=w)},
which we aim to model or approximate by a parametric function of w.


## The implemented functions

Two functions for the semiparametric case and two functions for the nonparametric case are provided in this package. "spOR" and "npOR" will be of primary interest to most users. The former (spOR) assumes a parametric model for the conditional odds ratio, while the latter (npOR) does not assume anything about the functional form of the conditional odds ratio but uses a parametric working model as an approximation. Thus, npOR allows for robust and correct inference even when the true conditional odds ratio is very complex, while still allowing for interpretable estimates based on a user-specified parametric working model. On the other hand, spOR will give incorrect, and possibly misleading, inference when the parametric model is incorrect (does not contain the true OR function).

### For the semiparametric case:

We assume the partially-linear logistic-link model
logit(P(Y=1|A=a,W=w)) = a*b^T f(w) + h(w)
where b^T f(w) is a correct parametric model for the log conditional odds ratio logOR(w), and h(w), which is necessarily equal to logit(P(Y=1|A=0,W=w)), is left unspecified and estimated with machine-learning (The Highly Adaptive Lasso, R package: hal9001).

"spOR" provides semiparametric estimates and inference for a user-specified parametric functional form of the conditional odds ratio. The estimates and inference are only correct if the parametric form is correctly specified (i.e. captures the true conditional odds ratio). This method is semiparametric since no further assumptions are made on any other feature of the data-generating probability distribution. Notably, the nuisance function P(Y=1|A=0,W), which is variation independent of the conditional odds ratio, is left totally unspecified (nonparametric) and is learned from the data using machine-learning. This function allows for outcome missingness, however for correct inference, the variables that predict both the outcome Y and outcome missingness Delta must come before the treatment A and be included in W. Rigorously, we require the the outcome Y and missingness mechanism Delta are conditionally independent conditional on A, W.

"spORMissing" is a more general semiparametric version of "spOR" that should be used if there are post-treatment variables Z that inform both the outcome and outcome missingness. For correct inference, this function requires the weaker assumption that  the outcome Y and missingness mechanism Delta are conditionally independent conditional on Z, A, W. Note this function still requires that the parametric form of the odds ratio is correctly specified for correct inference. 

### For the nonparametric case:

We make no assumptions on the functional form of the conditional odds ratio.

"npOR" provides nonparametric estimates and inference for a user-specified parametric working-model for the true conditional odds ratio. This is a nonparametric version of "spOR". Since the parametric form/model is only used to obtain a  best working-model approximation, this method provides correct estimates and inference without requiring the parametric model to be correct. Thus, this method is truly nonparametric. This function can handle the same type of missingness that "spOR" can.

"npORMissing" is a more general nonparametric version of "npOR" that should be used if there are post-treatment variables Z that inform both the outcome and outcome missingness. This is the nonparametric version of "spORMissing". Like "npOR", the parametric model is only used an approximation, so correct estimates and inference are obtainable without assuming the parametric model is correct. Thus, this method is truly nonparametric. This function can handle the same type of missingness that "spORMissing" can.


Each function allows the user to supply an arbitrary parametric form/working-model for the conditional log odds ratio. Estimates and 95% confidence intervals are returned. 

## Machine-learning and Targeted learning of the OddsRatio

Nuisance functions can be estimated using machine-learning (specifically the machine-learning pipeline R package tlverse/sl3), thereby avoiding biases due to model misspecification. To estimate the nuisance functions, for convenience, we allow the user to specify either an sl3 Learner object or a R formula object in which case the nuisance functions are parametrically estimated using glm (not recommended). By default, tlverse/hal9001 is used to estimate all nuisance functions. To allow for valid inference with black-box machine-learning, the initial estimators are debiased with the Targeted Maximum-Likelihood estimation (TMLE) framework, which gives compatible targeted substitution estimators that respect the constraints of the statistical model (unlike inferior estimating-equation-based methods like Double/Debiased Machine-Learning which do not respect any known information). TMLE is the state-of-the-art method for robust efficient semi/nonparametric inference, and even performs well at small sample sizes.


