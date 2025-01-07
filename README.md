# Balanced Boosting Algorithm for GAMLSS Using Adaptive Step Lengths

***

This repository provides R-functions for boosting generalized additive models for location, scale and shape (GAMLSS) 
with different step length methods as well as for conducting the corresponding cross-validation.
Currently only linear base-learners are included. 
The implemented step length methods are introduced in the paper 
*A Balanced Statistical Boosting Approach for GAMLSS via New Step Lengths* by Alexandra Daub, Andreas Mayr, 
Boyao Zhang and Elisabeth Bergherr. 

Available response variable distributions:
  * Gaussian
  * Negative Binomial
  * Weibull

Available step lengths methods:
  * fixed step lengths
  * numerically obtained shrunk optimal step lengths
  * analytically derived (approximated) shrunk optimal step lengths for one of the distribution parameters
  * base-learner based step lengths (update size is set to size of other update)

***

The repository contains **two main functions**:
  * boost_GAMLSS.R: boosting algorithm for fitting GAMLSS models with different step length methods
  * cv_risk.R: computation of the cross-validated risk (negative log-likelihood) for the boosted GAMLSS 

A **minimal working example** for the functions *boost_GAMLSS()* and *cv_risk()* is provided for the 
different response variable distributions:
  * example.R: 
      * simulate data
      * conduct a cross-validation
      * estimate model based on cross-validated stopping iteration
      * extract certain results from estimated model
