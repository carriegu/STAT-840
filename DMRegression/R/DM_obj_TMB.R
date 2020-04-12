#' Return the DMRegression TMB object constructed on the specific dataset
#'
#' @param Y Matrix of `n x J` responses.
#' @param X Matrix of `n x (P+1)` covariates.
#' @return The TMB object constructed on the specific dataset
#' @export
DM_obj_TMB <- function(Y, X, beta){

  # Extract some parameters
  n = nrow(Y)
  J = ncol(Y)
  P = ncol(X) - 1

  # Construct the TMB object
  tmb_model = MakeADFun(data=c(model = "DMRegressionFreq", #which model to use
                               list(X=X, Y=Y)), 
                        parameters=list(beta=beta),DLL="DMRegressionFreq_TMBExports")

  return(tmb_model)

}
