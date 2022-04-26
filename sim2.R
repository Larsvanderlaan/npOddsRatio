ns <- 2*c(   15000, 250, 500, 2500, 1000)
library(future)
library(future.apply)
library(MASS)
library(sl3)
library(R6)
true <- c(1)
plan(multiprocess)
fun <- function(n) {
  passes <- c()
  all_coefs<- c()
  for(i in 1:1000){





    X <- MASS::mvrnorm(n = n,
                       mu = c(0, 0 ),
                       Sigma = rbind(0.5*c(1, 0.5   ),
                                     c(0.5, 1 )
                       ))
    W1 <- bound(X[,1], c(-1.5,1.5))
    W2 <- bound(X[,2], c(-1.5,1.5))

    A <- rbinom(n, size = 1, plogis(0.75*(X %*% c(1,-1 ) )))
    Tvar <- rweibull(n,shape = 3, scale = 1/exp(0.1 * (W1 + W2 - 0.5)))
    Tcenter <- Tvar-1

    Q <- plogis( A * (1 +  Tcenter)  + W1 + W2 + Tcenter)

    J <- rbinom(n, 1, Q )

    R <- rbinom(n, size = 1, prob = plogis(0.5*(W1 + W2 + A - 1 + Tcenter)))
    R <- 1
    data <- data.frame(W1, W2, A, Tcenter, J, R)
    data <- data[data$R==1,]



    A <- data$A
    Y <- data$J
    W <- as.matrix(data[,c("Tcenter", "W1", "W2")])
    # library(future)
    #plan(multisession)

    task_A <- sl3_Task$new(data, covariates = c("Tcenter", "W1", "W2"), outcome = "A")
    task_Y <- sl3_Task$new(data, covariates = c("Tcenter", "W1", "W2", "A"), outcome = "J")
    data1 <- data
    data1$A <- 1
    data0 <- data
    data0$A <- 0
    task_Y0 <- sl3_Task$new(data0, covariates = c("Tcenter", "W1", "W2", "A"), outcome = "J")
    task_Y1 <- sl3_Task$new(data1, covariates = c("Tcenter", "W1", "W2", "A"), outcome = "J")

    lrnr_sp <- Lrnr_glm_semiparametric$new(~ 1 + Tcenter, lrnr_baseline = Lrnr_glmnet$new(), append_interaction_matrix = TRUE, family = binomial())
    lrnr_sp <- lrnr_sp$train(task_Y)

    EY1 <- lrnr_sp$predict(task_Y1)
    EY0 <- lrnr_sp$predict(task_Y0)

    lrnr_A <- Lrnr_glmnet$new()


    lrnr_A <- lrnr_A$train( task_A)

    pA1 <- pmin(pmax(lrnr_A$predict(task_A), 0.005), 1-0.005)


    # spout <- spglm( formula = ~1 + Tcenter, W = c("Tcenter", "W1", "W2"), A  = "A", Y = "J", data = data, estimand = "OR", sl3_Learner_A = lrnr_A, sl3_Learner_Y = lrnr_Y, append_interaction_matrix = TRUE )
    spout <- spOR(formula = ~1 + Tcenter, W = c("Tcenter", "W1", "W2"), A  = "A", Y = "J", data = data, EY1, EY0, pA1)

    # doMC::registerDoMC(cores = 11)
    # spout <- spOR(formula = ~1 + Tcenter, W, A, Y, data = data, Delta = NULL, sl3_learner_A = Lrnr_hal9001$new(max_degree = 2, num_knots = c(10,8), fit_control = list(parallel = TRUE)), smoothness_order_Y0W = 1, max_degree_Y0W = 2, num_knots_Y0W = c(10, 8),fit_control = list(parallel = TRUE))
    pout <- glm(J~ W1 + W2 + A * (1 + Tcenter) + Tcenter , family = binomial, data = data)


    coefs<-spout$coefs
    print(coefs[,1])

    lower <- coefs[,4]
    upper <- coefs[,5]
    print("ci width sp")
    print(upper - lower)
    pass <- lower <= true & upper >= true

    print(coef( pout)[c(4,6)])
    beta <-coef( pout)[c(4,6)]
    ci <- confint( pout)[c(4,6),,drop = F]
    lower <- ci[,1]
    upper <- ci[,2]
    print("ci width parametric")
    print(upper - lower)


    pass <- c(pass, lower <= true & upper >= true)
    all_coefs <- cbind(all_coefs, c(coefs,beta))
    passes <- cbind(passes, pass)
    print(rowMeans(passes))


  }
  return(list(n=n, all_coefs = all_coefs, passes = passes, true = true))


}

output <- future_lapply(ns, fun)
save(output, file =  "spORsim2.RData")


