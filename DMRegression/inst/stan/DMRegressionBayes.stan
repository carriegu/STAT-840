// @brief Dirichlet-Multinomial Regression with Stan.
//
// @details The DM model used here
//
// ```
// y(x)|Phi ~ Multinomial(n,Phi)
// Phi|Gamma ~ Dirichlet(Gamma)
// Gamma = exp(X %*% beta)
// beta follows flat prior

functions{
  // funtion to calcualte the DM - negative loglikelihood
  //' @param beta Matrix of `(P+1) x J` parameters.
  //' @param Y Matrix of `n x J` responses.
  //' @param X Matrix of `n x (P+1)` covariates.
  //' @return Scalar value of the negative loglikelihood.
  real negloglikelihood_lpdf(real[ ] Y, real [ ] X, int N, int P, int J, matrix beta){
    //define additional local variables
    real ll; //variable to hold log likelihood
    vector [J] Gamma_i; //gamma for observation i
    real Gamma_plus; // row_sum of gamma_i
    real n; // row_sum of Y_i 
    //for the variables appearing in sum loops, set them to zero
    
    for (j in 1:J){
      Gamma_i[j]=exp(to_row_vector(X)*beta[,j]);
    }

    Gamma_plus = sum(Gamma_i);
    n = sum(Y);


    ll=0;
    ll = ll + lgamma(n+1) + lgamma(Gamma_plus) - lgamma(n+Gamma_plus);
    for(j in 1:J){
      ll = ll + lgamma(Y[j] + Gamma_i[j]) - lgamma(Gamma_i[j]) - lgamma(Y[j]+1);
      }
    
    return (ll);
  }
}


// The input data is a matrix 'y' with `n x J` responses
data {
  int N; ///< Number of observations 
  int J; ///< Number of microbial taxon 
  int P; /// < Number of covariate 
  real X [N,(P+1)]; ///< `Matrix of `n x (P+1)` covariates.
  real Y[N,J]; ///< Matrix of `n x J` responses.
}

// The parameters accepted by the model. Our model
// accepts one parameter 'beta'.
parameters {
  matrix[(P+1),J] beta; ///< Matrix of `(P+1) x J` parameters.
}



// The model to be estimated. We model the output
// 'y' to be a DM distribution with parameters Gamma = exp(X%*%beta)
model {
  for (i in 1:N){
    Y[i]~negloglikelihood_lpdf(X[i,],N, P, J, beta);
  }

}

