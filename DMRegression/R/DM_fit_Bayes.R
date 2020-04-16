#' Estimate beta coefficients of the Dirichlet-Multinomial regression model
#' using Bayes approach with a flat prior.
#'
#' @param Y Matrix of `n x J` responses.
#' @param X Matrix of `n x (P+1)` covariates.
#' @param iterations Number of iterations for MCMC posterior sampeling.
#' @param chains Number of chains for MCMC posterior sampeling.
#' @return point estimation of `(P+1) x J` beta coefficients.
#' @return Samples via MCMC.
#' @export
DM_fit_Bayes = function(Y, X, iterations, chains){

  # Extract some parameters
  n = nrow(Y)
  J = ncol(Y)
  P = ncol(X) - 1

  # instantiate the model corresponding to p(mu, sigma, lambda | y, X) in Stan
  DM_data <- list(N=n,J=J,P=P, X = X, Y = Y)
  DM_fit <- rstan::sampling(stanmodels$DMRegressionBayes, data = DM_data, iter = iterations,
                     verbose = TRUE, chains = chains)

  DM_point_est <- rstan::get_posterior_mean(DM_fit) # The point estimation using posterior mean
  total_par_num = (P+1)*J
  DM_point_est <- c(DM_point_est)[c(1:total_par_num)]
  DM_point_est <- matrix(DM_point_est, nrow = P+1, ncol = J, byrow = TRUE)
  DM_samples <- rstan::extract(DM_fit) # The MCMC samples

  return(list(DM_point_est = DM_point_est, DM_samples = DM_samples))
}

