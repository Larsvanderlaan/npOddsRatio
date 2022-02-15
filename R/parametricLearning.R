

autoML <- function(n, p, fast_analysis = F) {
  if(p >= 25) {
    max_degree <- 1
  } else {
    max_degree <- 2
  }
  if(n<=50) {
    lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 1, num_knots =  c(1,0))
  }
  else if(n<=100) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 2, num_knots =  c(1,1))
    } else {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 1, num_knots =  c(1,1))
    }
  } else if(n<=250) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 2, num_knots =  c(7,3))
    } else {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(3,1))
    }
  } else if(n<=500) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 2, num_knots =  c(12,5))
    } else {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(10,5))
    }
  } else if(n<=1000) {
    lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(25,10))
  } else {
    lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(25,15))
  }
  return(lrnr)
}

causalGLM <- function(formula, W, A, Y,  learning_method = c("autoML", "glm", "glmnet", "gam", "xgboost"), method = c("CATE", "OR", "RR"), inference_type = c(  "semiparametric", "parametric"), weights = NULL, sl3_Learner = NULL, smoothness_order = 1, max_degree = 2, num_knots = c(20,10) ){
  method <- match.arg(method)
  inference_type <- match.arg(inference_type)
  W <- as.matrix(W)
  n <- nrow(W)
  p <- ncol(p)
  if(is.null(sl3_Learner)) {
    if(learning_method == "autoML" ) {
      sl3_Learner <- autoML(n,p)
    } else if(learning_method == "glmnet" ) {
      sl3_Learner <- Lrnr_glmnet$new()
    } else if(learning_method == "glm" ) {
      sl3_Learner <- Lrnr_glm$new()
    } else if(learning_method == "gam" ) {
      sl3_Learner <- Lrnr_gam$new()
    } else if(learning_method == "xgboost" ) {
      sl3_Learner <- Stack$new( Lrnr_glmnet$new(), Lrnr_xgboost$new(max_depth =3), Lrnr_xgboost$new(max_depth =4), Lrnr_xgboost$new(max_depth =5), Lrnr_xgboost$new(max_depth =6), Lrnr_xgboost$new(max_depth =7))
      sl3_Learner <- make_learner(Pipeline, Lrnr_cv$new(sl3_Learner), Lrnr_cv_selector$new(loss_squared_error))
    }
  }
  if(inference_type == "semiparametric") {
    if(method == "CATE") {
      return(spCATE(formula_CATE =  formula, W, A, Y,  sl3_Lrnr_A = sl3_Learner, sl3_Lrnr_Y = sl3_Learner,   weights = weights,  smoothness_order_Y0W = smoothness_order, max_degree_Y0W = max_degree, num_knots_Y0W = num_knots))
    }
    if(method == "OR") {
      return(spOR(formula_CATE =  formula, W, A, Y,  sl3_Lrnr_A = sl3_Learner, sl3_Lrnr_Y = sl3_Learner,   weights = weights,  smoothness_order_Y0W = smoothness_order, max_degree_Y0W = max_degree, num_knots_Y0W = num_knots))
    }

  }

}



