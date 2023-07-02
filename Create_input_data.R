#Modified from Michail Spitieris
create_data <-function(angle, 
                       t_theta, 
                       nc, 
                       n_theta, 
                       n_E, 
                       sigma_theta, 
                       sigma_E, 
                       prior_sigma_theta_mu, 
                       prior_sigma_theta_var, 
                       range_sigma_theta_upper, 
                       prior_sigma_delta_mu, 
                       prior_sigma_delta_var, 
                       prior_l_delta_mu, 
                       prior_l_delta_var, 
                       prior_R_mu, prior_R_var,
                       seed = 0, n_pred=0){
  t_E = seq(0, max(t_theta)-0.05, length.out = n_E)
  ytheta_real = angle
  yE_real = rep(0, length = n_E)
  
  #Add noise
  set.seed(seed)
  noise_theta = rnorm(n_theta*nc, 0, sigma_theta)
  noise_E = rnorm(n_E*nc, 0, sigma_E)
  ytheta_real = rep(ytheta_real,each=nc)
  yE_real = rep(yE_real,each=nc)
  
  
  ytheta_noisy = ytheta_real + noise_theta
  yE_noisy = yE_real + noise_E
  y_noisy = c(ytheta_noisy, yE_noisy)
  ytheta = ytheta_noisy
  yE = yE_noisy
  
  
  y = c(ytheta, yE)
  
  t_theta_pred = 0
  t_E_pred = 0
  
  
  data_noisy_pred = list(n_theta = nc*n_theta, n_E = nc*n_E, t_theta = matrix(rep(t_theta,each=nc),ncol=1)
                         , t_E = matrix(rep(t_E,each=nc),ncol=1), ytheta = ytheta, yE = yE
                         , prior_sigma_theta_mu = prior_sigma_theta_mu, prior_sigma_theta_var = prior_sigma_theta_var
                         , range_sigma_theta_upper =range_sigma_theta_upper, prior_sigma_delta_mu=prior_sigma_delta_mu, prior_sigma_delta_var=prior_sigma_delta_var
                         , prior_l_delta_mu=prior_l_delta_mu, prior_l_delta_var=prior_l_delta_var
                         , n_theta_pred = n_pred, n_E_pred=n_pred, t_theta_pred=matrix(rep(t_theta_pred, each=nc), ncol = 1), t_E_pred=matrix(rep(t_E_pred, each=nc), ncol = 1), prior_R_mu = prior_R_mu, prior_R_var=prior_R_var
  )
  
  # data_mod_true = data.frame(t_theta=t_theta, t_E=t_E, E = 0, Angle = angle)
  return(list(data_noisy_pred = data_noisy_pred, y=y))
  
}


# For predicting
create_data_pred <-function(angle, 
                            t_theta,
                            nc, 
                            n_theta, 
                            n_E, 
                            n_pred,
                            sigma_theta, 
                            sigma_E, 
                            prior_sigma_theta_mu, 
                            prior_sigma_theta_var, 
                            range_sigma_theta_upper, 
                            prior_sigma_delta_mu, 
                            prior_sigma_delta_var, 
                            prior_l_delta_mu, 
                            prior_l_delta_var, 
                            seed = 0){
  
  t_E = seq((t_theta[1]+t_theta[2])/2, (t_theta[n_theta]+t_theta[n_theta-1])/2, length.out = n_E)
  ytheta_real = angle
  yE_real = rep(0, length = length(t_E))
  
  
  #Add noise
  set.seed(seed)
  sigma_theta = rnorm(n_theta*nc, 0, sigma_theta)
  sigma_E = rnorm(n_E*nc, 0, sigma_E)
  ytheta_real = rep(ytheta_real,each=nc)
  yE_real = rep(yE_real,each=nc)
  
  ytheta_noisy = ytheta_real + sigma_theta
  yE_noisy = yE_real + sigma_E
  y_noisy = c(ytheta_noisy, yE_noisy)
  
  #If not normalizing
  ytheta = ytheta_noisy
  yE = yE_noisy
  
  
  y = c(ytheta, yE)
  
  t_theta_pred = seq(0.1, max(t_theta)-0.2, length=n_pred)
  num_t_preds = 1
  while(num_t_preds<=n_pred){
    prop = runif(1, 0, max(t_theta))
    print(prop)
    if (prop %in% t_theta){
    }
    else{
      t_theta_pred[num_t_preds] = prop
      num_t_preds = num_t_preds+1
    }
  }
  t_E_pred = seq(0.01, max(t_theta)-0.03, length=n_pred)
  
  for (i in 1:n_pred){
    for (j in 1:n_theta){
      if (abs(t_theta_pred[i]-t_theta[j])<=1e-2){
        print("Wopsi_theta")
        t_theta_pred[i] = t_theta_pred[i]+0.05
      }
      # else{
      #   print(t_theta_pred[i]-t_theta[j])
      # }
    }
  }
  for (i in 1:n_pred){
    for (j in 1:n_E){
      if (abs(t_theta_pred[i]-t_E[j])<=1e-2){
        print("Wopsi_E")
        t_theta_pred[i] = t_theta_pred[i]+0.05
      }
      # else{
      #   print(t_theta_pred[i]-t_theta[j])
      # }
    }
  }
  
  
  for (i in 1:n_pred){
    for (j in 1:n_theta){
      if (abs(t_theta_pred[i]-t_theta[j])<=1e-2){
        print("Wopsi")
        # t_theta_pred[i] = t_theta_pred[i]+0.01
      }
      # else{
      #   print(t_theta_pred[i]-t_theta[j])
      # }
    }
  }
  
  
  
  data_noisy_pred = list(n_theta = nc*n_theta, n_E = nc*n_E, t_theta = matrix(rep(t_theta,each=nc),ncol=1)
                         , t_E = matrix(rep(t_E,each=nc),ncol=1), ytheta = ytheta, yE = yE
                         , prior_sigma_theta_mu = prior_sigma_theta_mu, prior_sigma_theta_var = prior_sigma_theta_var
                         , range_sigma_theta_upper =range_sigma_theta_upper, prior_sigma_delta_mu=prior_sigma_delta_mu, prior_sigma_delta_var=prior_sigma_delta_var
                         , prior_l_delta_mu=prior_l_delta_mu, prior_l_delta_var=prior_l_delta_var
                         , n_theta_pred = n_pred, n_E_pred=n_pred, t_theta_pred=matrix(rep(t_theta_pred, each=nc), ncol = 1), t_E_pred=matrix(rep(t_E_pred, each=nc), ncol = 1), prior_R_mu = prior_R_mu, prior_R_var=prior_R_var
  )
  
  # data_mod_true = data.frame(t_theta=t_theta, t_E=t_E, I = 0, P = angle)
  return(list(data_noisy_pred = data_noisy_pred, y=y))
  
}

# For predicting
create_data_extrapolate <-function(angle, 
                                   t_theta,
                                   nc, 
                                   n_theta, 
                                   n_E, 
                                   n_pred,
                                   sigma_theta, 
                                   sigma_E, 
                                   prior_sigma_theta_mu, 
                                   prior_sigma_theta_var, 
                                   range_sigma_theta_upper, 
                                   prior_sigma_delta_mu, 
                                   prior_sigma_delta_var, 
                                   prior_l_delta_mu, 
                                   prior_l_delta_var, 
                                   seed = 0){
  
  t_E = seq((t_theta[1]+t_theta[2])/2, (t_theta[n_theta]+t_theta[n_theta-1])/2, length.out = n_E)
  ytheta_real = angle
  yE_real = rep(0, length = length(t_E))
  
  
  #Add noise
  set.seed(seed)
  sigma_theta = rnorm(n_theta*nc, 0, sigma_theta)
  sigma_E = rnorm(n_E*nc, 0, sigma_E)
  ytheta_real = rep(ytheta_real,each=nc)
  yE_real = rep(yE_real,each=nc)
  
  ytheta_noisy = ytheta_real + sigma_theta
  yE_noisy = yE_real + sigma_E
  y_noisy = c(ytheta_noisy, yE_noisy)
  
  #If not normalizing
  ytheta = ytheta_noisy
  yE = yE_noisy
  
  #If normalizing
  # tranform y in [0,1]
  # ytheta_temp = ytheta_noisy
  # yE_temp = yE_noisy
  # rthetaE = range(c(ytheta_temp, yE_temp))
  # ytheta = (ytheta_temp - rthetaE[1])/diff(rthetaE)
  # yE = (yE_temp - rthetaE[1])/diff(rthetaE)
  
  
  y = c(ytheta, yE)
  
  t_theta_pred = seq(max(t_theta)-0.2, max(t_theta)+3, length=n_pred)

  t_E_pred = seq(0.01, max(t_theta)-0.03, length=n_pred)
  # 

  # 
  
  
  data_noisy_pred = list(n_theta = nc*n_theta, n_E = nc*n_E, t_theta = matrix(rep(t_theta,each=nc),ncol=1)
                         , t_E = matrix(rep(t_E,each=nc),ncol=1), ytheta = ytheta, yE = yE
                         , prior_sigma_theta_mu = prior_sigma_theta_mu, prior_sigma_theta_var = prior_sigma_theta_var
                         , range_sigma_theta_upper =range_sigma_theta_upper, prior_sigma_delta_mu=prior_sigma_delta_mu, prior_sigma_delta_var=prior_sigma_delta_var
                         , prior_l_delta_mu=prior_l_delta_mu, prior_l_delta_var=prior_l_delta_var
                         , n_theta_pred = n_pred, n_E_pred=n_pred, t_theta_pred=matrix(rep(t_theta_pred, each=nc), ncol = 1), t_E_pred=matrix(rep(t_E_pred, each=nc), ncol = 1), prior_R_mu = prior_R_mu, prior_R_var=prior_R_var
  )
  
  # data_mod_true = data.frame(t_theta=t_theta, t_E=t_E, I = 0, P = angle)
  return(list(data_noisy_pred = data_noisy_pred, y=y))
  
}




# transform to the real scale
transform_post = function(y, fit){
  post_df = as.data.frame(extract(fit))
  return(post_df)
}

