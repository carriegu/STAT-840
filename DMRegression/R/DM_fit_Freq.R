#' Estimate beta coefficients of the Dirichlet-Multinomial regression model
#' using Frequentist approach.
#'
#' @param Y Matrix of `n x J` responses.
#' @param X Matrix of `n x (P+1)` covariates.
#' @return beta_hat the fitted matrix of `(P+1) x J` parameters and the corresponding Fisher-info.
#' @export
DM_fit_Freq <- function(Y, X){

  # Extract some parameters
  n = nrow(Y)
  J = ncol(Y)
  P = ncol(X) - 1

  #initial parameters
  beta0<-matrix(rnorm(n=(P+1)*J), nrow = P+1, ncol = J)
  #beta0<-matrix(0, nrow = P+1, ncol = J)

  # Construct the TMB object
  tmb_model = MakeADFun(data=c(model = "DMRegressionFreq", #which model to use
                               list(X=X, Y=Y)), 
                        parameters=list(beta=beta0),DLL="DMRegressionFreq2_TMBExports")


  #Optimization
  mle <- optim(par = c(beta0),
               fn = tmb_model$fn,
               gr = tmb_model$gr,
               method = "L-BFGS-B",
               control = list(maxit = 2000000,  ndeps = 1e-6),
               upper = 3, lower = -3)

  # The point estimation
  DM_point_est <- mle$par # The point estimation using mle
  total_par_num = (P+1)*J
  DM_point_est <- c(DM_point_est)[c(1:total_par_num)]
  DM_point_est <- matrix(DM_point_est, nrow = P+1, ncol = J)

  # The Fisher-Information Matrix
  Fisher_Info = tmb_model$he(c(mle$par))
  return(list(DM_point_est = DM_point_est, Fisher_Info = Fisher_Info))

}
