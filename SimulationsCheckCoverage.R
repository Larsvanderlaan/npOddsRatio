


library(simcausal)

k <- 0 # iterations
##### TEST COVERAGE FOR npORMissing

passes <- c()
for(i in 1:200){
  #n <- 100000
  # lrnr <- Lrnr_hal9001_custom$new(max_degree = 3, num_knots = c(30,20,1), fit_control = list(parallel = TRUE))
  # task <- sl3_Task$new(data, covariates =c("W1", "W2") ,outcome = "Y")
  # lrnr_Y0 <- lrnr$train(task[data$A==0])
  # Q0 <- lrnr_Y0$predict(task)
  # true <- coef(glm.fit(cbind(rep(1,n), data$W1), data$Y, offset = qlogis(Q0), weights = data$A, family  = binomial()))
  #

  n <- 3500
  D <- DAG.empty()

  D <- D +
    node("W1", distr = "runif", min = -1, max = 1) +
    node("W2", distr = "runif", min = -1, max = 1) +
    node("A", distr = "rbinom", size = 1, prob = plogis(W1 + W2) )+
    node("Z1", distr = "rgamma", shape=3, rate = 2) +
    node("Z", distr = "rconst", const=min(Z1,4)) +
    node("Q", distr = "rconst", const = plogis(-1+A + A*(-1.5*W1  -2*(Z-1.5)*(Z<=1.5) + 2*W1*Z   ) + 0.5*W1+ (Z-1.3) + (Z-1.5)*(Z>=1.5) *(1.5+W1) + (Z-2.5)*(Z>=2.5)+ 0.5*W2))+

    node("Y", distr = "rbinom", size = 1, prob = Q) +
    node("G", distr = "rconst",const = 1- 0.6*plogis(0.4*(W1 + W2  +A-0.5) + (Z-1.3)*(2.5+W1) + 2*(Z-2)*(Z>=2) ))+
    node("Delta", distr = "rbinom", size = 1, prob = G)

  setD <- set.DAG(D )
  data <- sim(setD, n = n)
  true <- c(1.3511903, 0.9620972)
  library(doMC)
  registerDoMC()
  lrnr <- Lrnr_hal9001_custom$new(max_degree = 3, num_knots = c(20,10,1), fit_control = list(parallel = TRUE))

  fit <- npORMissing(~ poly(W1,degree=1, raw = T)   , W = data[,c("W1", "W2")], A = data$A, Y = data$Y , Z = data$Z, sl3_learner_A = Lrnr_gam$new(),
                     Delta = data$Delta, sl3_learner_default = lrnr )
  ci <- fit$coefs

  passes <- cbind(passes, ci[,3] <= true & ci[,4] >= true)
  print(rowMeans(passes))

  print(apply(passes,1,table))
}



##### TEST COVERAGE FOR npOR

passes <- c()
for(i in 1:k){
  n <- 2500
  D <- DAG.empty()

  D <- D +
    node("W1", distr = "runif", min = -1, max = 1) +
    node("W2", distr = "runif", min = -1, max = 1) +
    node("A", distr = "rbinom", size = 1, prob = plogis(W1 + W2) )+
    node("Q", distr = "rconst", const = plogis(-1+A + A*(W1  ) + W1 + W2))+
    node("Y", distr = "rbinom", size = 1, prob = Q) +
    node("G", distr = "rconst",const = 1- 0.7*plogis(W1 + W2  +A-0.5 ))+
    node("Delta", distr = "rbinom", size = 1, prob = G)

  setD <- set.DAG(D )
  data <- sim(setD, n = n)
  true <- c(1,1)
  lrnr <- Lrnr_hal9001_custom$new(max_degree = 2, num_knots = c(10,1))


  fit <- npOR(~ poly(W1,degree=1, raw = T)   , W = data[,c("W1", "W2")], A = data$A, Y = data$Y ,
                     Delta = data$Delta, sl3_learner_default = lrnr )
  ci <- fit$coefs

  passes <- cbind(passes, ci[,3] <= true & ci[,4] >= true)
  print(rowMeans(passes))

  print(apply(passes,1,table))
}




##### TEST COVERAGE FOR spORMissing




k <- 1000
passes1 <- c()
passes2 <- c()
for(i in 1:k){
  n <- 3500
  D <- DAG.empty()
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "runif", min = -1, max = 1) +
    node("W2", distr = "runif", min = -1, max = 1) +
    node("A", distr = "rbinom", size = 1, prob = plogis(W1 + W2) )+
    node("Z1", distr = "rgamma", shape=3, rate = 2) +
    node("Z", distr = "rconst", const=min(Z1,4)) +
    node("Q", distr = "rconst", const = plogis(-1+A + A*(-1.5*W1  -2*(Z-1.5)*(Z<=1.5) + 2*W1*Z   ) + 0.5*W1+ (Z-1.3) + (Z-1.5)*(Z>=1.5) *(1.5+W1) + (Z-2.5)*(Z>=2.5)+ 0.5*W2))+

    node("Y", distr = "rbinom", size = 1, prob = Q) +
    node("G", distr = "rconst",const = 1- 0.7*plogis(0.4*(W1 + W2  +A-0.5) + (Z-1.3)*(2.5+W1) + 2*(Z-2)*(Z>=2) ))+
    node("Delta", distr = "rbinom", size = 1, prob = G)

  setD <- set.DAG(D )
  data <- sim(setD, n = n)
  data$Y <- data$Y * data$Delta
  true <- c(1.3104008, 0.9977911) # True coefficients estimated using very large sample and machine-learning algorithm
  library(doMC)
  registerDoMC(8)
  #fit <- npOR(~1 + W1^2 + W2^2, W = data[,c("W1", "W2")], A = data$A, Y = data$Y, glm_formula_A = ~ .^2, glm_formula_Y = ~ .^2  )
  #### INCORRECT as does not adjust for Z
  fit <- spOR(~ poly(W1,degree=1, raw = T)   , W = data[,c("W1", "W2")], A = data$A, Y = data$Y ,
                     Delta = data$Delta, num_knots_Y0W = c(30,10,1), max_degree_Y0W = 3,   sl3_learner_A = Lrnr_hal9001_custom$new(max_degree=2,num_knots = c(20,1)),fit_control = list(parallel = T)  )
  ci <- fit$coefs

  passes1 <- cbind(passes1, ci[,3] <= true & ci[,4] >= true)
  #### CORRECT as does not adjust for Z
  fit <- spORmissing(~ poly(W1,degree=1, raw = T)   , W = data[,c("W1", "W2")], A = data$A, Y = data$Y , Z =  data$Z,
                     Delta = data$Delta, num_knots_Y0W = c(30,10,1), max_degree_Y0W = 3, sl3_learner_Delta = Lrnr_hal9001_custom$new(max_degree=2,num_knots = c(20,1)), sl3_learner_A = Lrnr_hal9001_custom$new(max_degree=2,num_knots = c(20,1)),fit_control = list(parallel = T)  )
  ci <- fit$coefs

  passes2 <- cbind(passes2, ci[,3] <= true & ci[,4] >= true)

  # COMPARE
  print(rowMeans(passes1))
  print(rowMeans(passes2))
  print(apply(passes1,1,table))
  print(apply(passes2,1,table))
}




##### TEST COVERAGE FOR spOR





passes <- c()

for(i in 1:k){
  n <- 2500
  D <- DAG.empty()

  D <- D +
    node("W1", distr = "runif", min = -1, max = 1) +
    node("W2", distr = "runif", min = -1, max = 1) +
    node("A", distr = "rbinom", size = 1, prob = plogis(W1 + W2) )+
    node("Q", distr = "rconst", const = plogis(-1+A + A*(W1  ) + W1 + W2))+
    node("Y", distr = "rbinom", size = 1, prob = Q) +
    node("G", distr = "rconst",const = 1- 0.1*plogis(W1 + W2  +A-0.5 ))+
    node("Delta", distr = "rbinom", size = 1, prob = G)

  setD <- set.DAG(D )
  data <- sim(setD, n = n)
  data$Y <- data$Y * data$Delta
  true <- c(1, 1)
  library(doMC)
  registerDoMC(8)
  #fit <- npOR(~1 + W1^2 + W2^2, W = data[,c("W1", "W2")], A = data$A, Y = data$Y, glm_formula_A = ~ .^2, glm_formula_Y = ~ .^2  )
  #### INCORRECT as does not adjust for Z
  fit <- spOR(~ poly(W1,degree=1, raw = T)   , W = data[,c("W1", "W2")], A = data$A, Y = data$Y ,
                     Delta = data$Delta, num_knots_Y0W = c(30,10,1), max_degree_Y0W = 3, sl3_learner_A = Lrnr_hal9001_custom$new(max_degree=2,num_knots = c(20,1)),fit_control = list(parallel = T)  )

  ci <- fit$coefs

  passes <- cbind(passes, ci[,3] <= true & ci[,4] >= true)

  # COMPARE
  print(rowMeans(passes))

  print(apply(passes,1,table))

}




n <- 2500
D <- DAG.empty()

D <- D +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("W2", distr = "runif", min = -1, max = 1) +
  node("A", distr = "rbinom", size = 1, prob = plogis(W1 + W2) )+
  node("Q", distr = "rconst", const = plogis(-1+A + A*(W1  ) + W1 + W2))+
  node("Y", distr = "rbinom", size = 1, prob = Q) +
  node("G", distr = "rconst",const = 1- 0.1*plogis(W1 + W2  +A-0.5 ))+
  node("Delta", distr = "rbinom", size = 1, prob = G)

setD <- set.DAG(D )
data <- sim(setD, n = n)
data$Y <- data$Y * data$Delta
true <- c(1, 1)
library(doMC)
registerDoMC(8)
#fit <- npOR(~1 + W1^2 + W2^2, W = data[,c("W1", "W2")], A = data$A, Y = data$Y, glm_formula_A = ~ .^2, glm_formula_Y = ~ .^2  )
#### INCORRECT as does not adjust for Z
fit <- spOR(~ poly(W1, degree = 2, raw = T) + poly(W2, degree = 2, raw = T)    , W = data[,c("W1", "W2")], A = data$A, Y = data$Y ,
            Delta = data$Delta, num_knots_Y0W = c(30,10,1), max_degree_Y0W = 3, sl3_learner_A = Lrnr_hal9001_custom$new(max_degree=2,num_knots = c(20,1)),fit_control = list(parallel = T)  )

ci <- fit$coefs

summary(fit)


