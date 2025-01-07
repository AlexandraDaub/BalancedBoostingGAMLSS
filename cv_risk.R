library(parallel)

#' Cross-Validated Risk Calculation for Boosted GAMLSS Models
#'
#' This function computes the cross-validated risk (negative log-likelihood) for a 
#' Generalized Additive Model for Location, Scale, and Shape (GAMLSS) obtained via boosting.
#'
#' @param y A numeric vector of the response variable.
#' @param data A data frame or matrix of covariates. Each column corresponds to a covariate.
#' @param distribution A character string specifying the distribution to be used. 
#'        Options are "Gaussian" (default), "NegBinom" (Negative Binomial), or "Weibull".
#' @param m_stop An integer specifying the maximum number of boosting iterations. Default is 1000.
#' @param center_x A logical value. If `TRUE` (default), the predictors are centered (mean-subtracted).
#' @param folds A matrix of cross-validation folds. Each column represents a fold. 
#'        Defaults to 10-fold cross-validation.
#' @param methods_sl A character vector of length 2 specifying the step length methods for 
#'        the location (first) and scale (second) parameter. 
#'        Options for first distribution parameter: "LS" (Optimal, obtained numerically via line search, default), 
#'        "A" (Optimal, obtained based on analytical formula), "F" (Fixed), "BL" (Base-Learner based).
#'        Options for second distribution parameter: "LS" (default), "BL", "F".
#'        Note: "BL" has to be combined with "LS", "A" or "F".
#' @param searchInt_max A numeric vector of length 2 specifying the maximum search interval 
#'        for optimizing step lengths for location (first) and scale (second) parameter. Default is `c(10, 10)`.
#' @param c A numeric value specifying the fixed step length (used if methods_sl is "F"). Default is 1.
#' @param papply A function used for parallel processing. Default is `mclapply`.
#' @param show_progress A logical value. If `TRUE`, the starting time of each fold is printed. Default is `FALSE`.
#'
#' @return A list containing:
#'   \item{cvlike}{A matrix of cross-validated negative log-likelihoods for each fold (row) in each iteration (column).}
#'   \item{mstop}{An integer specifying the optimal number of boosting iterations based on the mean cross-validated log-likelihood.}
cv_risk <- function(y, data, distribution = "Gaussian", m_stop = 1000, center_x = T,
                    folds = cv(rep(1, length(y)), type = "kfold", B = 10),
                    methods_sl = c("LS", "LS"), searchInt_max = c(10, 10), c = 1,
                    papply = mclapply, show_progress = FALSE) {
  force(m_stop)
  force(folds)
  
  # perform cross-validation
  mod <- papply(1:ncol(folds), function(i) {
    if (show_progress){
      print(paste0("start fold ", i, " at: ", Sys.time()))
    }
    boost_GAMLSS(y, data, distribution = distribution, m_stop = m_stop, center_x = center_x, weights = as.logical(folds[, i]),
                 methods_sl = methods_sl, searchInt_max = c(searchInt_max[1], searchInt_max[2]), c = c)
  })
  
  # extract likelihoods
  likelihood <- -sapply(mod, `[[`, "like_test")
  
  output <- list()
  output$cvlike <- t(likelihood)
  output$mstop <- which.min(apply(output$cvlike, 2, mean, na.rm = T))
  return(output)
}

