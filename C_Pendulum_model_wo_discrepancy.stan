functions {
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
               real R) {
    int n_theta = rows(t_theta);
    int n_E = rows(t_E);
    matrix[n_theta + n_E, n_theta + n_E] K;

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
        K[i, n_theta + j] =-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) 
        + 0.100e1 * sigma * sigma * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) 
        + sigma * sigma * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) / R;

      }
    }
    
    
    // K_fu 
    
    for (i in 1:n_E){
      for (j in 1:n_theta){
       K[n_theta + i, j]  =-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) 
       + 0.100e1 * sigma * sigma * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) 
       + sigma * sigma * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) / R;

       // K[n_theta + i, j] = pow(sigma, 0.2e1) * K[n_theta + i, j];
     }
   }

   
   // K_ff
    for (i in 1:(n_E-1)){
     K[n_theta + i, n_theta +i] =0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) 
     - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.10000e1 * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) 
     - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.100e1 / R * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.100e1 * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)));

     for (j in (i+1):n_E){
      K[n_theta + i, n_theta +j] =0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) 
      - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[i,1] - t_E[j,1], 0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) 
      + 0.10000e1 * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) 
      - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) 
      + 0.100e1 / R * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) 
      + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) 
      + 0.100e1 * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) 
      + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)));

      // K[n_theta + i, n_theta +j] = pow(sigma, 0.2e1) * K[n_theta + i, n_theta +j];
      K[n_theta + j, n_theta +i] = K[n_theta + i, n_theta +j];
     }
     K[n_theta + n_E, n_theta +n_E] = 0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) 
     - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.10000e1 * sigma * sigma * pow(t_E[n_E,1] - t_E[n_E,1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) 
     - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.100e1 / R * sigma * sigma * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.100e1 * sigma * sigma * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)) 
     + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[n_E,1] - t_E[n_E,1], 0.2e1) * pow(l, -0.2e1)));

     }
     K[(n_theta + 1):(n_theta + n_E), (n_theta + 1):(n_theta + n_E)] = K[(n_theta + 1):(n_theta + n_E), (n_theta + 1):(n_theta + n_E)]
            + diag_matrix(rep_vector(pow(sigma_E, 0.2e1), n_E));
     return cholesky_decompose(K);
  }
} 


data {
  int<lower=1> n_theta;
  int<lower=1> n_E;
  // int<lower=1> n_theta_pred;
  // int<lower=1> n_E_pred;
  matrix[n_theta,1] t_theta;
  matrix[n_E,1] t_E;
  // matrix[n_theta_pred,1] t_theta_pred;
  // matrix[n_E_pred,1] t_E_pred;
  vector[n_theta] ytheta;
  vector[n_E] yE;
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
  // hyper-parameters
  real<lower=0.001> l;
  real<lower=0.001> sigma;
  // real mu_pendulum;
  real<lower=0, upper = range_sigma_theta_upper> sigma_theta;
  real<lower=0> sigma_E;
  // physical parameters
  real<lower=0> R;
}

model {
  // Chol. of PI kernel
  matrix[n_theta + n_E, n_theta + n_E] L_K = K_pendulum(t_theta, t_E, l, sigma, sigma_theta, sigma_E, R);
  // mean vector
  vector[n_theta + n_E] mu = mu_fn(n_theta, n_E); // Put in zero-vector + støy? sigma_E

  // Basis priors
  l ~ normal(1,1);
  sigma ~ normal(1,1);
  sigma_theta ~ normal(prior_sigma_theta_mu, prior_sigma_theta_var);
  sigma_E ~ normal(0, 0.1);

  // Basis priors
  R ~ normal(prior_R_mu, prior_R_var);



  y~ multi_normal_cholesky(mu, L_K);
}



