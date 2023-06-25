functions {
  
  int isPositiveDefinite(matrix M) {
    vector[rows(M)] eigenvalues = eigenvalues_sym(M);
    int numEigenvalues = num_elements(eigenvalues);
    
    for (i in 1:numEigenvalues) {
      if (eigenvalues[i] <= 0) {
        return 0;
      }
    }
    
    return 1;
  }

    vector mu_fn(int n_theta, 
               int n_E){
                 vector[n_theta]  m_theta = rep_vector(0,n_theta); 
                 vector[n_E]  m_E = rep_vector(0,n_E);
                 vector[n_theta + n_E] mu;
                 mu= append_row(m_theta, m_E);
                 return(mu);
  }
  matrix K_pendulum(matrix t_theta, 
               matrix t_E,
               real l,
               real sigma,
               real sigma_theta,
               real sigma_E,
               real R, 
               real sigma_delta,
               real l_delta) {
    int n_theta = rows(t_theta);
    int n_E = rows(t_E);
    matrix[n_theta + n_E, n_theta + n_E] K;
    matrix[n_E, n_E] K_d;
    

    // K_uu
    for (i in 1:(n_theta-1)){
      K[i,i] =   pow(sigma, 0.2e1);
      for (j in (i+1):n_theta){
        K[i,j] = exp(-0.5e0 * pow(t_theta[i,1] - t_theta[j,1], 0.2e1) * pow(l, -0.2e1));
        K[i,j] = pow(sigma, 0.2e1) * K[i,j];
        K[j,i] = K[i,j];
      }
      K[n_theta,n_theta] = pow(sigma, 0.2e1);
    }
    K[1:n_theta, 1:n_theta] = K[1:n_theta, 1:n_theta] + diag_matrix(rep_vector(pow(sigma_theta, 0.2e1), n_theta)); // observasjonsstøy
    
    
    // K_uf
    for (i in 1:n_theta){
      for (j in 1:n_E){
        K[i, n_theta + j] =-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) + sigma * sigma * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) / R;

      }
    }
    
    
    // K_fu 
    
    for (i in 1:n_E){
      for (j in 1:n_theta){
       K[n_theta + i, j]  =-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) + sigma * sigma * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) / R;

       // K[n_theta + i, j] = pow(sigma, 0.2e1) * K[n_theta + i, j];
     }
   }
   
     // Discrepancy K
    for (i in 1:(n_E-1)){
      K_d[i,i] = pow(sigma_delta, 0.2e1);
      for (j in (i+1):n_E){
        K_d[i,j] = exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l_delta, -0.2e1));
        K_d[i,j] = pow(sigma_delta, 0.2e1) * K_d[i,j];
        K_d[j,i] = K_d[i,j];
      }
      K_d[n_E,n_E] = pow(sigma_delta, 0.2e1);
    }

   
   // K_ff
    for (i in 1:(n_E-1)){
     K[n_theta + i, n_theta +i] =0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.10000e1 * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 / R * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)));

     for (j in (i+1):n_E){
      K[n_theta + i, n_theta +j] =0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[i,1] - t_E[j,1], 0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.10000e1 * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 / R * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)));

      // K[n_theta + i, n_theta +j] = pow(sigma, 0.2e1) * K[n_theta + i, n_theta +j];
      K[n_theta + j, n_theta +i] = K[n_theta + i, n_theta +j];
     }
     K[n_theta + n_E, n_theta +n_E] = 0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) + 0.10000e1 * sigma * sigma * pow(t_E[n_E,1] - t_E[n_E,1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 / R * sigma * sigma * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)));

     }
     K[(n_theta + 1):(n_theta + n_E), (n_theta + 1):(n_theta + n_E)] = K[(n_theta + 1):(n_theta + n_E), (n_theta + 1):(n_theta + n_E)]
            + diag_matrix(rep_vector(pow(sigma_E, 0.2e1), n_E)) + K_d[1:n_E, 1:n_E];
     return cholesky_decompose(K);
  }

matrix K_uu_pred(matrix t_theta_pred,
                 matrix t_theta,
                 real l,
                 real sigma
                 ){
    // Make K_uu for predictions. Takes in vector of observed angle times and angle times we wish to predict. Returns Kuu- matrix between those times with the estimated l and sigma
    int n_theta = rows(t_theta);
    int n_theta_pred = rows(t_theta_pred);
    matrix[n_theta_pred, n_theta] Kuu;

    for (i in 1:n_theta_pred){
      for (j in 1:n_theta){
        Kuu[i,j] = exp(-0.5e0*pow(t_theta_pred[i,1] - t_theta[j,1], 2) * pow(l, -2));
        Kuu[i,j] = pow(sigma, 2) * Kuu[i,j];
      }
    }
    return Kuu;

  }


  matrix K_uf_pred(matrix t_theta_pred,
                  matrix t_E,
                  real l,
                  real sigma,
                  real R
  ){
    // Make K_uf for predictions. Takes in vector of observed energy times and angle times we wish to predict. Returns Kuu- matrix between those times with the estimated R, l and sigma
   // Just changing the pred to not contain pred st the expression is written the same
    int n_theta = rows(t_theta_pred);
    matrix[n_theta,1] t_theta = t_theta_pred;
    int n_E = rows(t_E);
    matrix[n_theta, n_E] Kuf;

    // K_uf
    for (i in 1:n_theta){
      for (j in 1:n_E){
        Kuf[i, j] =-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1))
        + 0.100e1 * sigma * sigma * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1))
        + sigma * sigma * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) / R;
        }
      }
      return Kuf;
    }

vector Pgp_pred_rng(matrix t_theta,
                      matrix t_E,
                      vector ytheta,
                      vector yE,
                      matrix t_theta_pred,
                      real sigma,
                      real l,
                      real l_delta,
                      real sigma_delta,
                      real R,
                      real sigma_theta,
                      real sigma_E) {
    // Make mu_function and Covariance function from eqations in 2.6 predictions
    int n_theta = rows(t_theta);
    int n_E = rows(t_E);
    int N = n_theta + n_E;
    int n_theta_pred = rows(t_theta_pred);
    vector[n_theta_pred] f2;
    int is_pd;
    int is_pd2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[n_theta_pred, N] V_uT_pred;
      matrix[N, n_theta_pred] tri_L_K_times_V_uT_pred;
      vector[n_theta_pred] f2_mu;
      matrix[n_theta_pred, n_theta_pred] cov_f2;
      matrix[n_theta_pred, n_theta_pred] diag_delta;
      vector[N] y;
      y[1:n_theta] = ytheta[1:n_theta];
      y[(n_theta+1):N] = yE[1:n_E];
      
      
      // Make cholesky-decomposition of K(t,t') based on the observed data (with neasurenebt noise)
      L_K =  K_pendulum(t_theta, t_E, l, sigma, sigma_theta, sigma_E, R, sigma_delta, l_delta);
      // The left division of y by a lower-triangular view of K; algebraically equivalent to the less efficient and stable form inverse(tri(K)) * y, where tri(L_K) is the lower-triangular portion of L_K with the above-diagonal entries set to zero.
      K_div_y = mdivide_left_tri_low(L_K, y);
      
      // The right division of K_div_y' by a triangular view of L_K; algebraically equivalent to the less efficient and stable form K_div_y' * inverse(tri(L_K)), where tri(L_K) is the lower-triangular portion of L_K with the above-diagonal entries set to zero.
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';
      
      
      // This is V_u^*T containing [Kuu(t_u*,t_u), K_uf(t_u*, t_f)]
      V_uT_pred[1:n_theta_pred, 1:n_theta] = K_uu_pred(t_theta_pred, t_theta, l, sigma);
      V_uT_pred[1:n_theta_pred, (n_theta+1):N] = K_uf_pred(t_theta_pred, t_E, l, sigma, R);
      // Assuming mu = 0. From equation, we have V_u^*T K^(-1)y
      f2_mu = (V_uT_pred * K_div_y);
      // The left division of qT_p' by a triangular view of L_K; algebraically equivalent to the less efficient and stable form inverse(tri(L_K)) * qT_p', where tri(L_K) is the lower-triangular portion of L_K with the above-diagonal entries set to zero.
      tri_L_K_times_V_uT_pred = mdivide_left_tri_low(L_K, V_uT_pred');
      
      cov_f2 = K_uu_pred(t_theta_pred, t_theta_pred, l, sigma) - tri_L_K_times_V_uT_pred' * tri_L_K_times_V_uT_pred;

      diag_delta = diag_matrix(rep_vector(1e-4, n_theta_pred));
      
      // is_pd = isPositiveDefinite(tri_L_K_times_V_uT_pred' * tri_L_K_times_V_uT_pred);
      is_pd = isPositiveDefinite(cov_f2);
      f2 = multi_normal_rng(f2_mu, cov_f2+diag_delta);

    }
    return f2;
  }
}

 
data {
  int<lower=1> n_theta;
  int<lower=1> n_E;
  matrix[n_theta,1] t_theta;
  matrix[n_E,1] t_E;
  vector[n_theta] ytheta;
  vector[n_E] yE;
  
  int<lower=1> n_theta_pred;
  int<lower=1> n_E_pred;
  matrix[n_theta_pred,1] t_theta_pred;
  // matrix[n_E_pred,1] t_E_pred;
  // var is actually the value of sd, meaning that the variance is prior_sigma_delta_var^2
  real<lower=0>prior_sigma_delta_mu;
  real<lower=0>prior_sigma_delta_var;
  real<lower=0>prior_l_delta_mu;
  real<lower=0>prior_l_delta_var;
  real <lower=0> prior_sigma_theta_mu;
  real<lower=0>prior_sigma_theta_var;
  real<lower=0>range_sigma_theta_upper;
  real<lower=0>prior_R_mu;
  real<lower=0>prior_R_var;
}

transformed data {
  vector[n_theta + n_E] y = append_row(ytheta, yE);
}


parameters {
  // hythetaer-parameters
  real<lower=0.001> l;
  real<lower=0.001> sigma;
  // real mu_pendulum;
  real<lower=0, upper= range_sigma_theta_upper> sigma_theta;
  real<lower=0.001>sigma_delta;
  real<lower=0.001>l_delta;
  real<lower=0> sigma_E;
  // physical parameters
  real<lower=0> R;
  // real<lower=0, upper=1> C;
}



model {
  // Chol. of K
  matrix[n_theta + n_E, n_theta + n_E] L_K = K_pendulum(t_theta, t_E, l, sigma, sigma_theta, sigma_E, R, sigma_delta,l_delta);
  // mean vector
  // vector[n_theta + n_E] mu = rep_vector(0, n_theta+n_E); // Put in zero-vector + støy? sigma_E
  vector[n_theta + n_E] mu = mu_fn(n_theta, n_E); // Put in zero-vector + støy? sigma_E
  // Basis priors
  l ~ normal(1,1);
  sigma ~ normal(1,1);
  sigma_theta ~ normal(prior_sigma_theta_mu, prior_sigma_theta_var);
  l_delta ~ normal(prior_l_delta_mu,prior_l_delta_var);
  sigma_delta~ normal(prior_sigma_delta_mu,prior_sigma_delta_var);
  sigma_E ~ normal(0, 0.1);

  R ~ normal(prior_R_mu, prior_R_var);

  
  y~ multi_normal_cholesky(mu, L_K);
}



generated quantities {
  vector[n_theta_pred] f_theta;
  vector[n_theta_pred] y_theta;
  // Find the expected values
  f_theta = Pgp_pred_rng(t_theta, t_E, ytheta, yE, t_theta_pred, sigma, l, l_delta, sigma_delta, R, sigma_theta, sigma_E);
  // Add estimated measurement noise to the expected values for the measurements
  for (n1 in 1:n_theta_pred){
    y_theta[n1] = normal_rng(f_theta[n1], sigma_theta);
  }
}





