#' Calculate the negative loglikelihood of the Dirichlet-Multinomial regression model.
#'
#' @param beta Matrix of `(P+1) x J` parameters.
#' @param Y Matrix of `n x J` responses.
#' @param X Matrix of `n x (P+1)` covariates.
#' @return Scalar value of the negative loglikelihood.
#' @export
DM_negloglikelihood = function(beta, X, Y){

  ll = 0

  # Extract some parameters
  n = nrow(Y)
  J = ncol(Y)
  P = ncol(X) - 1

  # Reshape beta to a matrix
  beta = matrix(beta, nrow = P+1, ncol = J)

  Gamma = exp(X %*% beta)
  n = apply(Y, 1, sum)
  Gamma_plus = apply(Gamma, 1, sum)
  for(i in c(1:nrow(X))){
    ll = ll + lgamma(n[i]+1) + lgamma(Gamma_plus[i]) - lgamma(n[i]+Gamma_plus[i])
    for(j in c(1:ncol(Y))){
      ll = ll + lgamma(Y[i,j] + Gamma[i,j]) - lgamma(Gamma[i,j]) - lgamma(Y[i,j]+1)
    }
  }
  return(-ll)
}
