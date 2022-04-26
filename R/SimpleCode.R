
library(R6)
spOR <- function(formula, W, A, Y, data, EY1, EY0, pA1) {
  n <- length(A)
  data <- as.data.frame(data)
  A <- data[, A]
  Y <- data[,Y]
  V <- model.matrix(formula, data = data[,W])
  EY <- ifelse(A==1, EY1, EY0)
  denom <- pmax((1-EY1)*(EY0), 1e-8)
  num <- EY1*(1-EY0)
  OR <- num/denom
  logOR <- log(pmax(OR, 1e-8))
  beta <-  (coef(glm(EY~AV-1, family = binomial(), offset = qlogis(EY0), data = list(AV=A*V, A = A, EY = EY ))))

  EY1 <- plogis(qlogis(EY0) + V %*% beta)

  EY <- ifelse(A==1, EY1, EY0)

  # Compute initial estimate of variance of coefficient estimates
  h_star <-  as.vector(-(pA1 * EY1 * (1-EY1)) / (pA1*EY1*(1-EY1) + (1-pA1)*EY0*(1-EY0)))
  H_star <- V*(A + h_star)
  scale <- apply(V,2, function(v){colMeans_safe(as.vector(  EY1*(1-EY1) * EY0*(1-EY0) * pA1 * (1-pA1) / (pA1 * EY1*(1-EY1) + (1-pA1) *EY0*(1-EY0) )) * v*V)})
  scale_inv <- solve(scale)
  var_unscaled <- as.matrix(var(H_star*(Y-EY)))
  var_scaled_init <-  scale_inv %*% var_unscaled  %*%  t(scale_inv)

  # Perform targeting step
  for(i in 1:200) {

    h_star <-  as.vector(-(pA1 * EY1 * (1-EY1)) / (pA1*EY1*(1-EY1) + (1-pA1)*EY0*(1-EY0)))
    H_star <- V*(A + h_star)
    H_star1 <- V*(1 + h_star)
    H_star0 <- V* h_star
    offset <- qlogis(EY)
    scale <- apply(V,2, function(v){colMeans_safe(as.vector(  EY1*(1-EY1) * EY0*(1-EY0) * pA1 * (1-pA1) / (pA1 * EY1*(1-EY1) + (1-pA1) *EY0*(1-EY0) )) * v*V)})

    scale_inv <- solve(scale)

    score <- max(abs(colMeans_safe( H_star%*%scale_inv*as.vector(Y-EY)) ))

    if(abs(score) <= 1/n || abs(score) <= min(0.5,0.5*min(sqrt(diag(var_scaled_init))))/sqrt(n)/log(n)){
      converged_flag <- TRUE
      break
    }

    clev_reduced <- H_star%*%scale_inv
    dir <- colMeans_safe(H_star%*%scale_inv*as.vector(Y-EY))
    dir <- dir/sqrt(mean(dir^2))
    risk <- function(epsilon) {

      Qeps <- plogis(offset + clev_reduced%*%dir * epsilon)
      loss <- -1  * ifelse(Y==1, log(Qeps), log(1-Qeps))
      return(mean(loss))


    }
    optim_fit <- optim(
      par = list(epsilon = 0.01), fn = risk,
      lower = 0, upper = 0.01,
      method = "Brent"
    )
    eps <-  (scale_inv %*%dir) * optim_fit$par


    #eps <- coef(glm(Y~X-1, family = binomial(), weights = weights, offset = offset, data = list(Y = Y, X = H_star)))

    EY0 <- as.vector(plogis(qlogis(EY0) +  H_star0 %*% eps ))
    EY1 <-  as.vector(plogis(qlogis(EY1) +  H_star1 %*% eps))
    EY <- ifelse(A==1, EY1, EY0)


  }
  print("SCORE: ")
  print(score)

  # Get coefficients
  denom <- pmax((1-EY1)*(EY0), 1e-8)
  num <- EY1*(1-EY0)
  OR <- num/denom
  logOR <- log(pmax(OR, 1e-8))
  estimates <- as.vector(coef(glm(logOR~V-1, family = gaussian(), data = list(V=V, logOR = logOR ))))



  var_scaled <- var_scaled_init

  se <- sqrt(diag(var_scaled))
  Zvalue <- abs(sqrt(n) * estimates/se)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)
  ci <- cbind(estimates,   se/sqrt(n), se ,estimates - 1.96*se/sqrt(n), estimates + 1.96*se/sqrt(n), Zvalue, pvalue)
  colnames(ci) <- c("coefs",   "se/sqrt(n)", "se","lower_CI", "upper_CI", "Z-value", "p-value")
  output <- list(coefs = ci, var_mat = var_scaled, score = score)
  class(output) <- c("spOR","npOddsRatio")

  return(output)



}



#' Scalable Highly Adaptive Lasso (HAL) semiparametric
#' @importFrom R6 R6Class
#' @importFrom origami folds2foldvec
#' @importFrom stats predict quasibinomial
#'

Lrnr_hal9001_semiparametric <- R6Class(
  classname = "Lrnr_hal9001_semiparametric",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(formula_sp, family = NULL, interaction_variable = "A", ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "weights", "ids"),
    .train = function(task) {
      args <- self$params


      formula_sp <- args$formula_sp
      trt <- args$interaction_variable
      A <- task$data[[trt]]
      W <- as.matrix(task$X)
      W <- W[, setdiff(task$nodes$covariates, trt), drop = F]

      V <- model.matrix(formula_sp, as.data.frame(W))

      outcome_type <- self$get_outcome_type(task)
      args$X <- as.matrix(W)
      args$X_unpenalized <- as.matrix(A * V)
      args$Y <- as.numeric(as.character(outcome_type$format(task$Y)))

      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family()
      }

      if (!any(grepl("fit_control", names(args)))) {
        args$fit_control <- list()
      }
      args$fit_control$foldid <- origami::folds2foldvec(task$folds)

      if (task$has_node("id")) {
        args$id <- task$id
      }

      if (task$has_node("weights")) {
        args$fit_control$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }

      # fit HAL, allowing glmnet-fitting arguments
      other_valid <- c(
        names(formals(glmnet::cv.glmnet)), names(formals(glmnet::glmnet))
      )

      fit_object <- sl3:::call_with_args(
        hal9001::fit_hal, args,
        other_valid = other_valid
      )

      return(fit_object)
    },

    .predict = function(task = NULL) {
      args <- self$params
      formula_sp <- args$formula_sp
      trt <- args$interaction_variable
      A <- task$data[[trt]]
      W <- as.matrix(task$X)
      W <- W[, setdiff(task$nodes$covariates, trt)]
      V <- model.matrix(formula_sp, as.data.frame(W))
      X_unpenalized <- as.matrix(A * V)
      predictions <- stats::predict(
        self$fit_object,
        new_data = as.matrix(W),
        new_X_unpenalized = X_unpenalized
      )
      if (!is.na(safe_dim(predictions)[2])) {
        p <- ncol(predictions)
        colnames(predictions) <- sprintf("lambda_%0.3e", self$params$lambda)
      }
      return(predictions)
    },

    .required_packages = c("hal9001", "glmnet")
  )
)

colMeans_safe <- function(X) {
  X <- as.matrix(X)
  if(ncol(X)==1){
    return(mean(as.vector(X)))
  }
  return(colMeans(X))
}
bound <- function(x, b){
  pmax(pmin(x,1-b),b)
}
