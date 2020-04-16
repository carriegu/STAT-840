#' Return the DMRegression Stan object constructed on the specific dataset
#'
#' @param Y Matrix of `n x J` responses.
#' @param X Matrix of `n x (P+1)` covariates.
#' @param ... Additional arguments to pass to [rstan::sampling()].
#' @return The Stan object constructed on the specific dataset
#' @export
DM_obj_Stan = function(Y, X, iterations, chains, ...){

  # Extract some parameters
  n = nrow(Y)
  J = ncol(Y)
  P = ncol(X) - 1

  # instantiate the model corresponding to p(mu, sigma, lambda | y, X) in Stan
  DM_data <- list(N=n,J=J,P=P, X = X, Y = Y)
  DM_fit <- rstan::sampling(stanmodels$DMRegressionBayes, data = DM_data, iter = iterations,
                            verbose = TRUE, chains = chains, ...)

  return(DM_fit)
}

