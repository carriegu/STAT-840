# STAT840 Final Project
# Unit tests for the **Bayes** model `Dirichlet_Multinomial`.
# This test file consists of 2 parts:
# 1. whether the likelihood provided by Stan is the same with calculated in R
# 2. whether the coefficients estimated are optimized (optimCheck)

require(optimCheck)
require(dirmult)
require(testthat)
require(rstan)

ntest <- 10 # number of tests
tol <- .05 # tolerance level

max_err <- rep(NA, ntest)
for(ii in 1:ntest) {

  # simulate random data
  J <- sample(2:10, 1)
  P <- sample(1:10, 1)
  n <- sample(50:100, 1)
  X <- matrix(rnorm(n*P), n, P)
  one <- rep(1, n)
  X <- cbind(one, X)
  N <- sample(c(1:1000), n)
  # simulate a initial beta matrix in order to simulate Y
  beta0 <- matrix(rnorm((P+1)*J), (P+1), J)
  beta_true = beta0
  # reparametrization for data generation
  Gamma = exp(X %*% beta0)
  Gamma_plus = apply(Gamma, 1, sum)
  theta = 1/(Gamma_plus + 1)
  pi = apply(Gamma, 2, function(x) {x / theta})
  # Simulate the responses Y from package 'dirmult'
  Y = simPop(J = 1, n = N[1], pi = pi[1,], theta = theta[1])$data
  for(jj in c(2:n)){
    Y = rbind(Y, simPop(J = 1, n = N[jj], pi = pi[jj,], theta = theta[jj])$data)
  }

  # Systematic unit testing Part 1 (Likelihood check).
  # Construct the Stan object for this specific dataset
  ##
  ## Stan_obj <- DM_obj_Stan(X=X, Y=Y, chains=1, iterations=2000)
  Stan_obj <- DM_obj_Stan(X=X, Y=Y, chains=1, iterations=1,
                          algorithm = "Fixed_param")
  ##
  lp_r=NULL
  lp_stan=NULL
  # check that difference between R and Stan nll is the same constant
  for(i in 1:ntest) {
    # randomly generate the parameter values
    beta_test =  matrix(rnorm((P+1)*J), (P+1), J)
    Pars <- list(beta=beta_test)
    # log-posterior calculation in R
    lp_r[i]<-DM_negloglikelihood(beta=beta_test,X=X,Y=Y)
    # log-posterior calculation in Stan
    upars <- unconstrain_pars(object = Stan_obj, pars = Pars)
    lp_stan[i]<--(log_prob(object = Stan_obj, upars = upars, adjust_transform = FALSE))
  }
  diff_range = abs(diff(range(abs(lp_r-lp_stan))))
  expect_lt(diff_range,tol)

  # Systematic unit testing Part 2 (Optim Check).
  # for each test check that `min(abs_err, rel_err) < tol`
  # for each element of the solution beta_hat.

  if(FALSE) {
    # fit the DM regression
    Bayes_fit <- DM_fit_Bayes(Y=Y, X=X, iterations = 2000, chains = 1)
    beta_hat = Bayes_fit$DM_point_est

    # projection plots, except we don't plot, instead saving the output
    # of each plot to an S3 object of type `optproj`.
    oproj <- optim_proj(xsol = c(beta_hat),
                        fun = function(beta)
                          DM_negloglikelihood(Y = Y, X = X, beta = beta),
                        maximize = FALSE, plot = FALSE)

    # `diff` calculates the abs and rel error between the
    # candidate solution `xsol` and the minimum in each projection plot.
    # see ?diff.optcheck for details.
    err <- abs(diff(oproj)) # abs and rel error
    max_err[ii] <- max(pmin(err[,"abs"], err[,"rel"]))
    expect_lt(max_err[ii],tol)
  }
}
