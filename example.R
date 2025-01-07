# example

source("boost_GAMLSS.R")
source("cv_risk.R")


# Gaussian GAMLSS ####

## simulate data
set.seed(1234)
n <- 500

for (i in 1:6) {
  assign(paste0("x", i), runif(n, -1, 1))
}
data1 <- as.matrix(as.data.frame(mget(paste0("x", 1:6))))
X1 <- cbind(1, data1)

beta1 <- c(0, 1, 2, 0.5, -1, 0, 0)
beta2 <- c(2, 0, 0, 0.2, 0.1, -0.1, -0.2)

mu1 <- X1 %*% beta1
sigma <- exp(X1 %*% beta2)

y1 <- rnorm(n, mu1, sigma)


## conduct cross-validation
folds <- mboost::cv(rep(1, length(y1)), type = "kfold", B = 10)
cv_gaussian <- cv_risk(y1, data1, distribution = "Gaussian", m_stop = 500, 
              folds = folds, center_x = TRUE, methods_sl = c("A", "LS"), 
              searchInt_max = c(100, 10))

## get stopping iteration
cv_gaussian$mstop


## fit model
res_gaussian <- boost_GAMLSS(y1, data1, distribution = "Gaussian", m_stop = cv_gaussian$mstop,
                             center_x = TRUE, methods_sl = c("A", "LS"), searchInt_max = c(10, 10),
                             extend_output = TRUE)


## investigate results

### coefficients at mstop
res_gaussian$b_mu
res_gaussian$b_si


### path plots without intercept
matplot(t(res_gaussian$mu_mat[2:7,]), type = "l",
        xlab = "iteration",
        ylab = "coefficient",
        main = "Path Plot Mu for Gaussian GAMLSS")
matplot(t(res_gaussian$si_mat[2:7,]), type = "l",
        xlab = "iteration",
        ylab = "coefficient",
        main = "Path Plot Sigma for Gaussian GAMLSS")


### get & plot step lengths
res_gaussian$v_mu
res_gaussian$v_si

plot(
  1:length(res_gaussian$v_mu), 
  res_gaussian$v_mu, 
  col = res_gaussian$v_mu_var,
  xlab = "iteration",
  ylab = "step length",
  main = "Step Lengths Mu for Gaussian GAMLSS")
legend("topright",
       legend = unique(res_gaussian$v_mu_var),
       col = unique(res_gaussian$v_mu_var),
       title = "Variable",
       pch = 1)

plot(
  1:length(res_gaussian$v_si), 
  res_gaussian$v_si, 
  col = res_gaussian$v_si_var,
  xlab = "iteration",
  ylab = "step length",
  main = "Step Lengths Sigma for Gaussian GAMLSS")
legend("topright",
       legend = unique(res_gaussian$v_si_var),
       col = unique(res_gaussian$v_si_var),
       title = "Variable",
       pch = 1)




# Negative Binomial GAMLSS ####

## simulate data
set.seed(1234)
n <- 500

for (i in 1:3) {
  assign(paste0("x", i), runif(n, -1, 1))
}
for(i in 4:6){
  assign(paste0("x", i), rbinom(n, 1, 0.5))
}
data2 <- as.matrix(as.data.frame(mget(paste0("x", 1:6))))
X2 <- cbind(1, data2)

beta1 <- c(-0.5, -0.5, 0.3, 0, 0.5, -0.3, 0)
beta2 <- c(0, 0, 0.6, -0.6, 0, -0.4, 0.4)

mu2 <- exp(X2 %*% beta1)
alpha <- exp(X2 %*% beta2)

y2 <- rnbinom(n, mu = mu2, size = 1/alpha)

## conduct cross-validation
folds <- mboost::cv(rep(1, length(y2)), type = "kfold", B = 10)
cv_negbinom <- cv_risk(y2, data2, distribution = "NegBinom", m_stop = 500, 
                       folds = folds, center_x = TRUE, methods_sl = c("A", "BL"), 
                       searchInt_max = c(20, 200))

## get stopping iteration
cv_negbinom$mstop

## fit model
res_negbinom <- boost_GAMLSS(y2, data2, distribution = "NegBinom", m_stop = cv_negbinom$mstop,
                             center_x = TRUE, methods_sl = c("A", "BL"), searchInt_max = c(20, 200))

## get coefficient estimates at mstop
res_negbinom$b_mu
res_negbinom$b_al




# Weibull GAMLSS ####

## simulate data
set.seed(1234)
n <- 500

for (i in 1:6) {
  assign(paste0("x", i), runif(n, -1, 1))
}
data3 <- as.matrix(as.data.frame(mget(paste0("x", 1:6))))
X3 <- cbind(1, data3)

beta1 <- c(0.6, 0.15, -0.2, 0.4, -0.25, 0, 0)
beta2 <- c(0, 0, 0, -0.15, 0.15, -0.1, 0.1)

lambda <- exp(X3 %*% beta1)
k <- exp(X3 %*% beta2)

y3 <- rweibull(n, scale = lambda, shape = k)

## conduct cross-validation
folds <- mboost::cv(rep(1, length(y3)), type = "kfold", B = 10)
cv_weibull <- cv_risk(y3, data3, distribution = "Weibull", m_stop = 700, 
                      folds = folds, center_x = TRUE, methods_sl = c("F", "F"))

## get stopping iteration
cv_weibull$mstop

## fit model
res_weibull <- boost_GAMLSS(y3, data3, distribution = "Weibull", m_stop = cv_weibull$mstop,
                            center_x = TRUE, methods_sl = c("F", "F"))

## get coefficient estimates at mstop
res_weibull$b_la
res_weibull$b_k
