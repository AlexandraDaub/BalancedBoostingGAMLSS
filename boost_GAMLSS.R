library(MASS)

#' Boosting a Generalized Additive Models for Location, Scale, and Shape (GAMLSS)
#'
#' This function implements a boosting algorithm for fitting GAMLSS models with different step length methods. 
#' Gaussian, Negative Binomial, and Weibull distributed response variables are supported.
#'
#' @param y A numeric vector of the response variable.
#' @param data A data frame or matrix of covariates. Each column corresponds to a covariate.
#' @param distribution A character string specifying the distribution to be used. 
#'        Options are "Gaussian" (default), "NegBinom" (Negative Binomial), or "Weibull".
#' @param m_stop An integer specifying the maximum number of boosting iterations. Default is 1000.
#' @param center_x A logical value. If `TRUE` (default), the predictors are centered (mean-subtracted).
#' @param weights An optional logical vector indicating a training and test data split. 
#'        If `NULL` (default), all data are used for both training and testing.
#' @param methods_sl A character vector of length 2 specifying the step length methods for 
#'        the location (first) and scale (second) parameter. 
#'        Options for first distribution parameter: "LS" (Optimal, obtained numerically via line search, default), 
#'        "A" (Optimal, obtained based on analytical formula), "F" (Fixed), "BL" (Base-Learner based).
#'        Options for second distribution parameter: "LS" (default), "BL", "F".
#'        Note: "BL" has to be combined with "LS", "A" or "F".
#' @param searchInt_max A numeric vector of length 2 specifying the maximum search interval
#'        for optimizing step lengths for location (first) and scale (second) parameter. Default is `c(10, 10)`.
#' @param c A numeric value specifying the fixed step length (used if methods_sl is "F"). Default is 1.
#' @param extend_output A logical value. If `TRUE`, additional outputs (e.g., coefficients estimates and 
#'        best fitted base-learners for both parameters in each iteration) are returned. Default is `FALSE`.
#'
#' @return A list containing:
#'   \item{b_param1, b_param2}{Estimated coefficients for the location (first) and scale (second) parameter.}
#'   \item{eta_param1, eta_param2}{Predictors for the location (first) and scale (second) parameter on the training set at mstop.}
#'   \item{v_param1, v_param2}{Step lengths for the location (first) and scale (second) parameter in each iteration.}
#'   \item{param1vparam2}{Distribution parameter updated in each iteration, 
#'         1 for location (first) parameter and 2 for scale (second) parameter.}
#'   \item{v_param1_var, v_param2_var}{Variable selected for update of the location (first) and scale (second) parameter in each iteration.}
#'   \item{like_train, like_test}{Log-likelihood values for the training and test set in each iteration.}
boost_GAMLSS <- function(y, data, distribution = "Gaussian", m_stop = 1000, center_x = TRUE, weights = NULL,
                         methods_sl = c("LS", "LS"), searchInt_max = c(10, 10), c = 1, extend_output = FALSE) {
  
  # check methods_sl input
  if (methods_sl[1]=="BL" & methods_sl[2]=="BL"){
    stop("Not all parameters can use the BL-ratio as step length\n")
  }
  
  if(methods_sl[2] == "A"){
    stop("No analytical optimal step lengths exist for second parameter\n")
  }
  
  # split data into training and test sets. If weights=NULL, then data = data_train = data_test
  if (!is.null(weights)) {
    weights <- as.logical(weights)
    y_train <- y[weights]
    data_train <- data[weights, ]
    y_test <- y[!weights]
    data_test <- data[!weights, ]
  }
  else {
    y_train <- y
    data_train <- data
    y_test <- y
    data_test <- data
  }
  
  if (center_x) {
    data_train <- apply(data_train, FUN = function(x) scale(x, scale = FALSE, center = TRUE), MARGIN = 2)
    data_test <- apply(data_test, FUN = function(x) scale(x, scale = FALSE, center = TRUE), MARGIN = 2)
  }
  else {
    data_train <- as.matrix(data_train)
    data_test <- as.matrix(data_test)
  }
  
  # add intercept
  X_train <- cbind(1, data_train)
  X_test <- cbind(1, data_test)
  p <- ncol(data_train)
  
  # define base-learners
  X_B <- lapply(1:p, function(j){
    cbind(1, data_train[, j])
  })
  
  # prepare fitting of base-learners
  B <- lapply(X_B, function(X_b){
    solve(t(X_b) %*% X_b) %*% t(X_b)
  })
  
  # initialize beta_param1 and beta_param2
  if(distribution == "Gaussian"){
    beta_param1 <- matrix(c(mean(y), rep(0, p)), ncol = 1)
    beta_param2 <- matrix(c(log(sd(y)), rep(0, p)), ncol = 1)
  } else if(distribution == "NegBinom"){
    beta_param1 <- matrix(c(log(mean(y)), rep(0, p)), ncol = 1)
    beta_param2 <- matrix(c(log((var(y)-mean(y))/mean(y)), rep(0, p)), ncol = 1)
  } else if(distribution == "Weibull"){
    empirical_distr_param <- fitdistr(y, densfun = "weibull")$estimate
    beta_param1 <- matrix(c(empirical_distr_param[2], rep(0, p)), ncol = 1)
    beta_param2 <- matrix(c(empirical_distr_param[1], rep(0, p)), ncol = 1)
  }
  
  # initialize other outputs
  param1vparam2 <- rep(0, m_stop)
  
  param1_train <- param1_test <- list()
  param2_train <- param2_test <- list()
  ngr_param1 <- list()
  ngr_param2 <- list()
  SSH_param1_train <- list()
  SSH_param2_train <- list()
  H_param1_train <- list()
  H_param2_train <- list()
  
  param1_mat <- param2_mat <- matrix(nrow = ncol(X_train), ncol = m_stop)
  
  v_param1 <- v_param1_var <- vector(length = m_stop)
  v_param2 <- v_param2_var <- vector(length = m_stop)
  
  eta_param1_train <- X_train %*% beta_param1
  eta_param1_test <- X_test %*% beta_param1
  eta_param2_train <- X_train %*% beta_param2
  eta_param2_test <- X_test %*% beta_param2
  step.length.param1 <- step.length.param2 <- 0
  
  like_train <- like_test <- vector(length = m_stop)
  
  # boosting start
  for (m in 1:m_stop) {
    if(methods_sl[1] == "BL"){ # start with update param2 if base-learner based step lengths used for param1 
      # potential update of parameter 2 ####
      # calculate negative gradient
      if(distribution == "Gaussian"){
        ngr_param2[[m]] <- (-1 + exp(-2 * eta_param2_train) * ((y_train - eta_param1_train)^2)) 
      } else if(distribution == "NegBinom"){
        ngr_param2[[m]] <- 1/exp(eta_param2_train)*log(1+exp(eta_param2_train)*exp(eta_param1_train)) + 
          (y_train-exp(eta_param1_train))/(1+exp(eta_param2_train)*exp(eta_param1_train)) - 
          digamma(y_train+1/exp(eta_param2_train))/exp(eta_param2_train) + digamma(1/exp(eta_param2_train))/exp(eta_param2_train)
      } else if(distribution == "Weibull"){
        ngr_param2[[m]] <- 1 + exp(eta_param2_train) * (log(y_train) - eta_param1_train) * (1 - exp( exp(eta_param2_train) * (log(y_train)-eta_param1_train) ) )
      }
      
      # fit base-learners
      B.coeff <- lapply(B, function(A){
        A %*% ngr_param2[[m]]
      })
      
      # compute SSR of each fit
      B.SSR <- lapply(1:length(B), function(j){
        sum((X_B[[j]] %*% B.coeff[[j]] - ngr_param2[[m]])^2)
      })
      
      # select best performing base-learner
      best_param2 <- which.min(B.SSR)
      v_param2_var[m] <- best_param2
      beta_best_param2 <- B.coeff[[best_param2]]
      h_param2_train <- X_train[,best_param2+1]*beta_best_param2[2] + beta_best_param2[1]
      h_param2_test <- X_test[,best_param2+1]*beta_best_param2[2] + beta_best_param2[1]
      SSH_param2_train[m] <- sum(h_param2_train^2)
      H_param2_train[[m]] <- h_param2_train
      
      # find step length
      if (methods_sl[2] == "LS") {
        phi_param2 <- function(v) {
          if(distribution == "Gaussian"){
            ret <- -sum(dnorm(y_train, eta_param1_train, exp(eta_param2_train + v * h_param2_train), log = TRUE))
          } else if(distribution == "NegBinom"){
            ret <- -sum(dnbinom(y_train, mu = exp(eta_param1_train), size = 1/exp(eta_param2_train + v * h_param2_train), log = TRUE))
          } else if(distribution == "Weibull"){
            ret <- -sum(dweibull(y_train, scale = exp(eta_param1_train), shape = exp(eta_param2_train + v * h_param2_train), log = TRUE))
          }
          return(ret)
        }
        v <- optimize(phi_param2, interval = c(-1, searchInt_max[2]))$min
      } else if(methods_sl[2] == "BL"){
        SSH_param1 <- sum(h_param1_train^2)
        SSH_param2 <- sum(h_param2_train^2)
        v <- v_param1[m] * SSH_param1/SSH_param2
      } else if (methods_sl[2] == "F") {
        v <- c
      }
      
      # the adaptive step length is 10% of the optimal step length
      v_param2[m] <- v
      step.length.param2 <- .1 * v_param2[m]
      
      # update coefficients
      beta_param2_st <- beta_param2
      beta_param2_st[1] <- beta_param2[1]+step.length.param2*beta_best_param2[1]
      beta_param2_st[1+best_param2] <- beta_param2[1+best_param2]+step.length.param2*beta_best_param2[2]
      
    }
    
    # potential update of parameter 1 ####
    # calculate negative gradient
    if(distribution == "Gaussian"){
      ngr_param1[[m]] <- (1/exp(eta_param2_train)^2) * (y_train - eta_param1_train)
    } else if(distribution == "NegBinom"){
      ngr_param1[[m]] <- (1/(1+exp(eta_param2_train)*exp(eta_param1_train))) * (y_train - exp(eta_param1_train))
    } else if(distribution == "Weibull"){
      ngr_param1[[m]] <- exp(eta_param2_train) * (exp( exp(eta_param2_train) * (log(y_train)-eta_param1_train) ) - 1)
    }
    
    # fit base-learners
    B.coeff <- lapply(B, function(A){
      A %*% ngr_param1[[m]]
    })
    
    # compute SSR of each fit
    B.SSR <- lapply(1:length(B), function(j){
      sum((X_B[[j]] %*% B.coeff[[j]] - ngr_param1[[m]])^2)
    })
    
    # select the best performing base-learner
    best_param1 <- which.min(B.SSR)
    v_param1_var[m] <- best_param1
    beta_best_param1 <- B.coeff[[best_param1]]
    h_param1_train <- X_train[,best_param1+1]*beta_best_param1[2] + beta_best_param1[1]
    h_param1_test <- X_test[,best_param1+1]*beta_best_param1[2] + beta_best_param1[1]
    SSH_param1_train[m] <- sum(h_param1_train^2)
    H_param1_train[[m]] <- h_param1_train
    
    # find step length
    phi_param1 <- function(v) {
      if(distribution == "Gaussian"){
        ret <- -sum(dnorm(y_train, eta_param1_train + v * h_param1_train, exp(eta_param2_train), log = TRUE))
      } else if(distribution == "NegBinom"){
        ret <- -sum(dnbinom(y_train, mu = exp(eta_param1_train + v * h_param1_train), size = 1/exp(eta_param2_train), log = TRUE))
      } else if(distribution == "Weibull"){
        ret <- -sum(dweibull(y_train, scale = exp(eta_param1_train + v * h_param1_train), shape = exp(eta_param2_train), log = TRUE))
      }
      return(ret)
    }
    if (methods_sl[1] == "LS") {
      v <- optimize(phi_param1, interval = c(-1, searchInt_max[1]))$min
    } else if (methods_sl[1] == "A") {
      if(distribution == "Gaussian"){
        b <- sum(h_param1_train^2)
        a <- sum(h_param1_train^2/exp(eta_param2_train)^2)
        v <- b / a
      } else if(distribution == "NegBinom"){
        if (m == 1){
          v <- optimize(phi_param1, interval = c(-1, searchInt_max[1]))$min                                                      
        } else {
          pos <- which(v_param1_var[1:m-1] == v_param1_var[m])
          if (length(pos) == 0){
            v_prev <- v_param1[[m-1]]
          } else {
            last_pos <- pos[length(pos)]
            v_prev <- v_param1[last_pos]
          }
          update_vprev <- exp(eta_param1_train + v_prev*h_param1_train)
          
          denom <- 1 + exp(eta_param2_train) * update_vprev
          
          a <- sum( h_param1_train^2 * update_vprev / denom )
          b <- sum( (h_param1_train*(y_train - (update_vprev * (1 - h_param1_train*v_prev) ))) / denom )
          
          v <- b/a
        }
      } else if(distribution == "Weibull"){
        if (m == 1){
          v <- optimize(phi_param1, interval = c(-1, searchInt_max[1]))$min
        } else {
          pos <- which(v_param1_var[1:m-1] == v_param1_var[m])
          if (length(pos) == 0){
            v_prev <- v_param1[[m-1]]
          } else {
            last_pos <- pos[length(pos)]
            v_prev <- v_param1[last_pos]
          }            
          update_vprev <- exp(-exp(eta_param2_train)*eta_param1_train - v_prev*exp(eta_param2_train)*h_param1_train)
          
          zeta <- exp(eta_param2_train) * h_param1_train * y_train^exp(eta_param2_train)
          
          a <- sum(-exp(eta_param2_train) * h_param1_train * zeta * update_vprev)
          b <- -sum(h_param1_train*exp(eta_param2_train)) + sum(zeta * update_vprev * (1+exp(eta_param2_train)*h_param1_train*v_prev))
          
          v <- -b/a
        }
      }
    } else if(methods_sl[1] == "BL"){
      SSH_param1 <- sum(h_param1_train^2)
      SSH_param2 <- sum(h_param2_train^2)
      v <- v_param2[m] * SSH_param2/SSH_param1
    } else if (methods_sl[1] == "F") {
      v <- c
    }
    
    # the adaptive step length is 10% of the optimal step length
    v_param1[m] <- v
    step.length.param1 <- .1 * v_param1[m]
    
    # update coefficients
    beta_param1_st <- beta_param1
    beta_param1_st[1] <- beta_param1[1]+step.length.param1*beta_best_param1[1]
    beta_param1_st[1+best_param1] <- beta_param1[1+best_param1]+step.length.param1*beta_best_param1[2]
    
    if(methods_sl[1] != "BL"){
      # potential update of parameter 2 ####
      # calculate negative gradient
      if(distribution == "Gaussian"){
        ngr_param2[[m]] <- (-1 + exp(-2 * eta_param2_train) * ((y_train - eta_param1_train)^2)) 
      } else if(distribution == "NegBinom"){
        ngr_param2[[m]] <- 1/exp(eta_param2_train)*log(1+exp(eta_param2_train)*exp(eta_param1_train)) + 
          (y_train-exp(eta_param1_train))/(1+exp(eta_param2_train)*exp(eta_param1_train)) - 
          digamma(y_train+1/exp(eta_param2_train))/exp(eta_param2_train) + digamma(1/exp(eta_param2_train))/exp(eta_param2_train)
      } else if(distribution == "Weibull"){
        ngr_param2[[m]] <- 1 + exp(eta_param2_train) * (log(y_train) - eta_param1_train) * (1 - exp( exp(eta_param2_train) * (log(y_train)-eta_param1_train) ) )
      }
      
      # fit base-learners
      B.coeff <- lapply(B, function(A){
        A %*% ngr_param2[[m]]
      })
      
      # compute SSR of each fit
      B.SSR <- lapply(1:length(B), function(j){
        sum((X_B[[j]] %*% B.coeff[[j]] - ngr_param2[[m]])^2)
      })
      
      # select the best performing base-learner
      best_param2 <- which.min(B.SSR)
      v_param2_var[m] <- best_param2
      beta_best_param2 <- B.coeff[[best_param2]]
      h_param2_train <- X_train[,best_param2+1]*beta_best_param2[2] + beta_best_param2[1]
      h_param2_test <- X_test[,best_param2+1]*beta_best_param2[2] + beta_best_param2[1]
      SSH_param2_train[m] <- sum(h_param2_train^2)
      H_param2_train[[m]] <- h_param2_train
      
      # find step length
      if (methods_sl[2] == "LS") {
        phi_param2 <- function(v) {
          if(distribution == "Gaussian"){
            ret <- -sum(dnorm(y_train, eta_param1_train, exp(eta_param2_train + v * h_param2_train), log = TRUE)) 
          } else if(distribution == "NegBinom"){
            ret <- -sum(dnbinom(y_train, mu = exp(eta_param1_train), size = 1/exp(eta_param2_train + v * h_param2_train), log = TRUE))
          } else if(distribution == "Weibull"){
            ret <- -sum(dweibull(y_train, scale = exp(eta_param1_train), shape = exp(eta_param2_train + v * h_param2_train), log = TRUE))
          }
          return(ret)
        }
        v <- optimize(phi_param2, interval = c(-1, searchInt_max[2]))$min
      } else if(methods_sl[2] == "BL"){
        SSH_param1 <- sum(h_param1_train^2)
        SSH_param2 <- sum(h_param2_train^2)
        v <- v_param1[m] * SSH_param1/SSH_param2
      } else if (methods_sl[2] == "F") {
        v <- c
      }
      
      # the adaptive step length is 10% of the optimal step length
      v_param2[m] <- v
      step.length.param2 <- .1 * v_param2[m]
      
      # update coefficients
      beta_param2_st <- beta_param2
      beta_param2_st[1] <- beta_param2[1]+step.length.param2*beta_best_param2[1]
      beta_param2_st[1+best_param2] <- beta_param2[1+best_param2]+step.length.param2*beta_best_param2[2]
    }
    
    # find and update the best best predictor ----
    # compare the positive likelihood
    if(distribution == "Gaussian"){
      like_param1 <- sum(dnorm(y_train, eta_param1_train + step.length.param1 * h_param1_train, exp(eta_param2_train), log = TRUE)) 
      like_param2 <- sum(dnorm(y_train, eta_param1_train, exp(eta_param2_train + step.length.param2 * h_param2_train), log = TRUE))
    } else if(distribution == "NegBinom"){
      like_param1 <- sum(dnbinom(y_train, mu = exp(eta_param1_train + step.length.param1 * h_param1_train), size = 1/exp(eta_param2_train), log = TRUE))
      like_param2 <- sum(dnbinom(y_train, mu = exp(eta_param1_train), size = 1/exp(eta_param2_train + step.length.param2 * h_param2_train), log = TRUE))
    } else if(distribution == "Weibull"){
      like_param1 <- sum(dweibull(y_train, scale = exp(eta_param1_train + step.length.param1 * h_param1_train), shape = exp(eta_param2_train), log = TRUE))
      like_param2 <- sum(dweibull(y_train, scale = exp(eta_param1_train), shape = exp(eta_param2_train + step.length.param2 * h_param2_train), log = TRUE))
    }
    
    # if like_param1 > like_param2, then update eta_param1
    if (like_param1 > like_param2) {
      beta_param1 <- beta_param1_st
      eta_param1_train <- eta_param1_train + step.length.param1 * h_param1_train 
      eta_param1_test <- eta_param1_test + step.length.param1 * h_param1_test
      param1vparam2[m] <- 1
    } else { # else update eta_param2
      beta_param2 <- beta_param2_st
      eta_param2_train <- eta_param2_train + step.length.param2 * h_param2_train 
      eta_param2_test <- eta_param2_test + step.length.param2 * h_param2_test
      param1vparam2[m] <- 2
    }
    
    # store updated coefficients, update and store distribution parameters and likelihood ----
    param1_mat[, m] <- beta_param1
    param2_mat[, m] <- beta_param2
    if(distribution == "Gaussian"){
      param1_train[[m]] <- eta_param1_train
      param1_test[[m]] <- eta_param1_test
      param2_train[[m]] <- exp(eta_param2_train)
      param2_test[[m]] <- exp(eta_param2_test)
      like_train[m] <- sum(dnorm(y_train, param1_train[[m]], param2_train[[m]], log = TRUE))
      like_test[m] <- sum(dnorm(y_test, param1_test[[m]], param2_test[[m]], log = TRUE))
    } else if(distribution == "NegBinom"){
      param1_train[[m]] <- exp(eta_param1_train)
      param1_test[[m]] <- exp(eta_param1_test)
      param2_train[[m]] <- exp(eta_param2_train)
      param2_test[[m]] <- exp(eta_param2_test)
      like_train[m] <- sum(dnbinom(y_train, mu = param1_train[[m]], size = 1/param2_train[[m]], log = TRUE))
      like_test[m] <- sum(dnbinom(y_test, mu = param1_test[[m]], size = 1/param2_test[[m]], log = TRUE))
    } else if(distribution == "Weibull"){
      param1_train[[m]] <- exp(eta_param1_train)
      param1_test[[m]] <- exp(eta_param1_test)
      param2_train[[m]] <- exp(eta_param2_train)
      param2_test[[m]] <- exp(eta_param2_test)
      like_train[m] <- sum(dweibull(y_train, scale = param1_train[[m]], shape = param2_train[[m]], log = TRUE))
      like_test[m] <- sum(dweibull(y_test, scale = param1_test[[m]], shape = param2_test[[m]], log = TRUE))
    }
  }
  
  # format the output
  rownames(beta_param1) <- c("(Intercept)", colnames(data))
  rownames(beta_param2) <- c("(Intercept)", colnames(data))
  
  if(distribution == "Gaussian"){
    name_param1 <- "mu"
    name_param2 <- "si"
  } else if(distribution == "NegBinom"){
    name_param1 <- "mu"
    name_param2 <- "al"
  } else if(distribution == "Weibull"){
    name_param1 <- "la"
    name_param2 <- "k" 
  }
  
  if(!extend_output){
    output <- list(b_param1 = t(beta_param1),
                   b_param2 = t(beta_param2),
                   eta_param1 = eta_param1_train,
                   eta_param2 = eta_param2_train,
                   v_param1 = v_param1,
                   v_param2 = v_param2,
                   param1vparam2 = param1vparam2,
                   v_param1_var = v_param1_var,
                   v_param2_var = v_param2_var,
                   like_train = like_train,
                   like_test = like_test,
                   call = sys.call())
    
    names(output) <- c(paste0("b_", name_param1),
                       paste0("b_", name_param2),
                       paste0("eta_", name_param1),
                       paste0("eta_", name_param2),
                       paste0("v_", name_param1),
                       paste0("v_", name_param2),
                       paste0(name_param1, "v", name_param2),
                       paste0("v_", name_param1, "_var"),
                       paste0("v_", name_param2, "_var"),
                       "like_train",
                       "like_test",
                       "call")
  } else{
    output <- list(b_param1 = t(beta_param1),
                   b_param2 = t(beta_param2),
                   param1_mat = param1_mat,
                   param2_mat = param2_mat,
                   eta_param1 = eta_param1_train,
                   eta_param2 = eta_param2_train,
                   SSH_param1_train = SSH_param1_train, 
                   SSH_param2_train = SSH_param2_train, 
                   H_param1_train = H_param1_train,
                   H_param2_train = H_param2_train,
                   v_param1 = v_param1,
                   v_param2 = v_param2,
                   param1vparam2 = param1vparam2,
                   v_param1_var = v_param1_var,
                   v_param2_var = v_param2_var,
                   like_train = like_train,
                   like_test = like_test,
                   call = sys.call())
    
    names(output) <- c(paste0("b_", name_param1),
                       paste0("b_", name_param2),
                       paste0(name_param1, "_mat"),
                       paste0(name_param2, "_mat"),
                       paste0("eta_", name_param1),
                       paste0("eta_", name_param2),
                       paste0("SSH_", name_param1, "_train"),
                       paste0("SSH_", name_param2, "_train"),
                       paste0("H_", name_param1, "_train"),
                       paste0("H_", name_param2, "_train"),
                       paste0("v_", name_param1),
                       paste0("v_", name_param2),
                       paste0(name_param1, "v", name_param2),
                       paste0("v_", name_param1, "_var"),
                       paste0("v_", name_param2, "_var"),
                       "like_train",
                       "like_test",
                       "call")
  }
  
  
  
  return(output)
}
