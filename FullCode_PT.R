# Load packages
library(foreign)
library(xtable)
library(stargazer)
library(rstan)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gmodels)
library(pracma)
library(BayesFactor)
library(MASS)
library(gridExtra)
library(bayesplot)
library(coda)
library(viridis)
library(shinystan)
library(philentropy)
library(FlexReg)
library(truncnorm)
# shinystan::launch_shinystan_demo()

# Set ggplot theme
theme_set(theme_classic())
# Define color palette

# my_colors <- brewer.pal(n = 2, name = "Dark2")


theme_set(
  theme_bw() + 
    theme(
      text = element_text(family = "Helvetica"),
      axis.title = element_text(size = rel(1.8)),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      plot.title = element_text(size = 25, vjust = 0, margin= margin(b = 30))
    )
)

# Set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = 3) # allocate 3 cores (for each model we run 3 chains in parallel)

#PENDULUM SIMULATION

# 1. LP: Linearized process
pendulum <- function(R, theta0, runtime, nummeasurements){
  times = seq(0.0, runtime,length.out = nummeasurements) 
  # times = runif(nummeasurements, 0,4)
  angle = theta0*cos(sqrt(1/R)*times)
  df = data.frame(Time = times, Angle = angle)
}

# 2. LPgp: Linearized process with added gaussian process
# gp: Generate Gaussian process values
gaussian_process <- function(l_delta, sigma_delta, runtime, seed) {
  
  set.seed(seed)
  # Define the kernel function
  kernel <- function(t, tp) {
    sigma_delta*exp(-0.5 * (t - tp)^2 / l_delta^2)
  }
  # Create the covariance matrix
  x <- seq(0, runtime, length.out = 100*runtime)
  cov_matrix <- matrix(0, nrow = 100*runtime, ncol = 100*runtime)
  for (i in 1:(100*runtime)) {
    for (j in 1:(100*runtime)) {
      cov_matrix[i, j] <- kernel(x[i], x[j])
    }
  }

  # Generate the Gaussian process values
  mu <- rep(0, (100*runtime))
  y <- mvrnorm(1, mu, cov_matrix)
  # Return a data frame with x and y columns
  data.frame(x = x, y = y)
}


# LPgp
simulate_pendulum_with_gp <- function(l_delta, sigma_delta, runtime, nummeasurements, R_real, theta0, seed){
  # Generate Gaussian process data

  gp <- gaussian_process(l_delta, sigma_delta, runtime, seed)
  #How many data points we want
  gpx <- gp$x[seq(1, length(gp$x), length.out = nummeasurements)]
  gpy <- gp$y[seq(1, length(gp$y), length.out = nummeasurements)]
  

  # Create a dataframe with the x and y data
  df <- data.frame(x = gpx, y = gpy)
  
  # Create the ggplot object
  plot <- ggplot(data = gp, aes(x, y)) +
    geom_line(color = "blue", linewidth = 1) +  # Adjust line color and thickness
    labs(title = "Gaussian Process", x = "t (s)", y = expression(delta(t)))
  
  # Save the plot as a PNG file
  ggsave("gp1.pdf")
  
  # Generate pendulum data
  pendulum_list <- pendulum(R_real, theta0, runtime, nummeasurements)
  

  # Add GP values to pendulum angle
  pendulum_list$Angle <- pendulum_list$Angle + gpy
  
  return(pendulum_list)
}


#3. TP: True idealized process
pendulum_true <- function(R_real,theta0, runtime, h, nummeasurements) {
  L=R_real*g
  t = c( seq(0, runtime, by=h))
  n = length(t)
  
  y=c(rep(0, n))
  v=c(rep(0, n))
  
  accel<-function(theta){
    # return(-(m*g*L/I)*sin(theta)) #-1/Rsin(theta)
    return(-(1/R_real)*sin(theta))
  }
  
  y[1] = theta0
  v[1] = 0
  
  for (i in 1:n){
    k1y = h*v[i]
    k1v = h*accel(y[i])
    
    k2y = h*(v[i]+0.5*k1v)
    k2v = h*accel(y[i]+0.5*k1y)
    
    k3y = h*(v[i]+0.5*k2v)
    k3v = h*accel(y[i]+0.5*k2y)
    
    k4y = h*(v[i]+k3v)
    k4v = h*accel(y[i]+k3y)
    
    # Update next value of y
    y[i+1] = y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0
    v[i+1] = v[i] + (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0
  }
  y=y[1:n]
  
  #How many data points we want
  t <- t[seq(1, length(t), length.out = nummeasurements)]
  y <- y[seq(1, length(y), length.out = nummeasurements)]
  
  df_true <- data.frame(Time=t, Angle=y)
  return(df_true)
}



#4. DP: damped process
pendulum_damped <- function(L, damping, g, theta0, omega0, nummeasurements, runtime, h) {
  # # Calculate the time step
  # dt <- runtime/(nummeasurements-1)
  
  
  # Create empty vectors to store the time and angle measurements
  times <- seq(0, runtime, by = h)
  angles <- numeric(length(times))
  
  # Implement the damped pendulum equation using the Euler method
  theta <- theta0
  omega <- omega0
  for (i in seq_along(times)) {
    angles[i] <- theta
    theta <- theta + omega * h
    omega <- omega - (damping*omega + (g/L)*sin(theta)) * h
  }
  
  times <- times[seq(1, length(times), length.out = nummeasurements)]
  angles <- angles[seq(1, length(angles), length.out = nummeasurements)]
  # Create a data frame with the time and angle measurements
  df <- data.frame(Time = times, Angle = angles)
  
  return(df)
}


#For energy calculation
calculate_energy <- function(df, R, m, g, L, theta0){
  L=R*g
  potential = m*g*L*(1-cos(df$Angle))
  kinetic =  m*g*L*(cos(df$Angle)-cos(theta0))
  Energy = potential + kinetic
  df = data.frame(Time=df$Time, Energy=Energy)
  return(df)
}

find_height <- function(df, R){
  L=R*g
  h=L*(1-cos(df$Angle))
  df = data.frame(Time = df$Time, Height = h)
}

#Function to plot pendulum
makeplot <- function(data, xaxis, yaxis){
  ggplot(data = data, aes(xaxis, yaxis)) + 
    theme_bw() + 
    geom_point(col = "tomato1", size = 2) + 
    labs(title = plottitle, subtitle = plotsubtitle) + 
    labs(x = "Time(s)", y = expression(paste("Angle ", theta, " (rad)"))) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5))
}


#Plot pendulum with and without noise
compare_with_noise <- function(sigma_theta){
  df=pendulum_true(R_real = 0.2, theta0 = pi/5, runtime = 4, h =0.0001, nummeasurements = 20)
  noisy = create_data(df$Angle, df$Time, nc = 1, n_theta = as.integer(length(df$Time)),
                      n_E = as.integer(length(df$Time)), sigma_theta, sigma_E, prior_sigma_theta_mu,
                      prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var,
                      prior_l_delta_mu, prior_l_delta_var, prior_R_mu, prior_R_var, seed) 
  noisydf=data.frame(Time=noisy$data_noisy_pred$t_theta, Angle=noisy$data_noisy_pred$ytheta)
  df=pendulum_true(R_real=0.2, theta0=pi/5, runtime=4, h=0.0001, 500)
  plot = ggplot(data = df, aes(x=Time,y = Angle,  color="True")) +
    geom_line(size=1) +
    scale_x_continuous(name = "Time (s)", limits = c(0, runtime)) +
    scale_y_continuous(name = "Angle (rad)") +
    geom_point(data = noisydf, aes(x= Time, y =Angle, color="Observations"), shape=4, stroke=2) +
    # labs(color = "Process") +
    ggtitle(expression(paste("Observations with ", sigma[u], "=0.1")))+
    scale_color_manual(name = " ",
                       values = c("Observations" = "black","True"="darkgreen"),
                       labels = c("Observations", "True"))
  ggsave("W_and_WO_noise01.pdf", plot, height = 5, width = 10)
}




#Fit model once
fit_model_once <- function(noisydf, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, 
                           prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var,model_type = "with_disc", seed) {
  
  if (model_type == "w_disc") {
    stanfile <- stanfile_w_disc
  } else if (model_type == "w_E_disc") {
    stanfile <- stanfile_w_E_disc
  } else if (model_type == "wo_disc"){
    stanfile <- stanfile_wo_disc
  } else if (model_type == "p_wo_disc"){
    stanfile = stanfile_P_wo_disc
  }
  
  result <- stan(file=stanfile,
                 data=noisydf$data_noisy_pred,
                 chains=3,
                 iter=numiters,
                 seed=seed)
  
  return(result)
}

#Save traceplot
save_trace <- function(result, name){
  # color_scheme_set("teal")
  color_scheme_set("mix-blue-red")
  theme_set(theme_classic())
  traceplot = mcmc_trace(result)
  traceplot
  ggsave(name, traceplot)
}


#Make plot of posterior histogram
create_posterior_histogram <- function(df, model, histtitle, variables, priors_mu = NULL, priors_var = NULL, plot_priors = TRUE){
  
  if (!is.null(priors_mu)) {
    priors_mu <- setNames(priors_mu, variables)
    priors_var <- setNames(priors_var, variables)
  }
  
  # Make histogram ready list
  y <- df$y # the original observed data
  pp_se <- transform_post(y, model)
  pl_df <- pp_se [, variables]
  pl_df$sample <- 1:nrow(pl_df)
  m_pl_df <- melt(pl_df, id="sample")
  
  # Create histogram and density plot
  histogram <- ggplot(data=m_pl_df) +
    ggtitle(histtitle) +
    geom_histogram(aes(x=value, y=..density.., color=variable), bins=50, inherit.aes = FALSE, fill="lavender") +
    #geom_density(aes(x=value, color=variable), size=1, alpha=0.4) +  # Add density plot
    theme_classic() +
    facet_wrap(~variable, nrow = 1, scales = "free") +
    theme(legend.position="none") +
    geom_vline(data=realval[realval$variable %in% variables,], aes(xintercept=real, color=variable), color = "springgreen4", size=1) +
    theme(text = element_text(size = 15))
  
  if (plot_priors & !is.null(priors_mu)) {
    for (var in variables) {
      # Get the prior mean and variance for the variable
      prior_mu <- priors_mu[var]
      prior_var <- priors_var[var]
      
      # Calculate the lower and upper bounds of the 75% confidence interval
      prior_lower <- qnorm(0.40, prior_mu, sqrt(prior_var))
      prior_upper <- qnorm(0.60, prior_mu, sqrt(prior_var))
      
      # Create a sequence of values within the 75% confidence interval
      seq_values <- seq(prior_lower, prior_upper, length.out = 100)
      
      # Create a data frame with the values and variable name
      prior_df <- data.frame(value = seq_values, variable = var) 
      
      # Add a column with the density values for the prior distribution
      prior_df$density <- dnorm(seq_values, mean = prior_mu, sd = sqrt(prior_var))* (seq_values > 0)
      
      # Add the prior distribution to the histogram
      histogram <- histogram + geom_line(data = prior_df, aes(x = value, y = density),
                                         color = "darkorchid4", size = 1, linetype = "dashed")
    }
  }
  
  return(histogram)
}

# Create combined histogram of posteriors
create_combined_histogram <- function(hist1, hist2, hist3) {
  grid.arrange(hist1, hist2, hist3, nrow = 3)
}

#Make tables of bias and credible interval width
create_table <- function(result_wo_disc, result_w_disc, result_w_E_disc, R_real) {
  
  # Define function to calculate bias and CI width
  calculate_bias_and_ci_width <- function(fit, R_real) {
    bias <- summary(fit)$summary["R", "mean"] - R_real
    ci_width <- summary(fit)$summary["R", "97.5%"] - summary(fit)$summary["R", "2.5%"]
    return(c(bias, ci_width))
  }
  
  # Create data frame with bias and CI width for each model
  df <- data.frame(
    Model = c("Without discrepancy", "Discrepancy on angle", "Discrepancy on energy"),
    Bias = c(calculate_bias_and_ci_width(result_wo_disc, R_real)[1], calculate_bias_and_ci_width(result_w_disc, R_real)[1], calculate_bias_and_ci_width(result_w_E_disc, R_real)[1]),
    CI_width = c(calculate_bias_and_ci_width(result_wo_disc, R_real)[2], calculate_bias_and_ci_width(result_w_disc, R_real)[2], calculate_bias_and_ci_width(result_w_E_disc, R_real)[2])
  )
  
  # Create LaTeX table
  latex_table <- xtable(df, digits = 4, caption = "Summary")
  
  # Print the LaTeX code for the table
  print(latex_table, include.rownames = FALSE, caption.placement = "top", caption.width = "0.9\\textwidth")
}

#Save the plots
save_plots <- function(combined_histogram, histogram_w_disc, histogram_wo_disc, histogram_w_E_disc, num) {
  ggsave(paste("combined_histogram", num, ".pdf"), plot = combined_histogram)
  ggsave(paste("w_disc_histogram", num, ".pdf"), plot = histogram_w_disc)
  ggsave(paste("wo_disc_histogram", num, ".pdf"), plot = histogram_wo_disc)
  ggsave(paste("w_E_disc_histogram", num, ".pdf"), plot = histogram_w_E_disc)
  # ggsave(paste("p_wo_disc_histogram", num, ".pdf"), plot = histogram_p_wo_disc)
}

#Run the fit several times, with same input data, but different seed for each computer model run
run_multiple_fits <- function(df, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, 
                              prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type, num_runs = 10) {
  
  # create an empty list to store the results
  results_list <- list()
  
  
  
  if (model_type == "w_disc") {
    stanfile <- stanfile_w_disc
  } else if (model_type == "w_E_disc") {
    stanfile <- stanfile_w_E_disc
  } else if (model_type == "wo_disc"){
    stanfile <- stanfile_wo_disc
  } else if (model_type == "p_wo_disc"){
    stanfile = stanfile_P_wo_disc
  }
  
  print(numiters)
  # loop over the desired number of runs
  for (i in 1:num_runs) {
    # set a new random seed for each run
    seed <- i*1000
    
    df2 <- create_data(df$Angle, df$Time, nc = nc, n_theta = as.integer(length(df$Time)),
                       n_E = as.integer(length(df$Time)), sigma_theta, sigma_E, prior_sigma_theta_mu,
                       prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var,
                       prior_l_delta_mu, prior_l_delta_var, prior_R_mu, prior_R_var, seed) 
    
    
    # Create plot of df2$data_noisy_pred$t_theta and df2$data_noisy_pred$ytheta using ggplot
    p_theta <-   ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_theta, y = ytheta)) +
      geom_point() +
      xlab("Time") +
      ylab("Theta")
    plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$ytheta)
    print("yo")
    # Save plot as pdf file
    ggsave("noisy_prediction_theta.pdf", p_theta)
    
    # Create plot of df2$data_noisy_pred$t_E and df2$data_noisy_pred$yE using ggplot
    p_E <-  ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_E, y = yE)) +
      geom_point() +
      xlab("Time") +
      ylab("E")
    
    # Save plot as pdf file
    ggsave("noisy_prediction_E.pdf", p_E)
    # run the model
    result <- stan(file=stanfile,
                   data=df2$data_noisy_pred,
                   chains=3,
                   iter=numiters,
                   seed=seed)
    fit_summary <- summary(result)
    print(names(fit_summary))
    print(fit_summary)
    print(fit_summary$summary)
    # mu_tau_summary <- summary(fit_summary, pars = c("R", "sigma_theta"), probs = c(0.01, 0.99))$summary
    # print(mu_tau_summary)
    
    name = paste(model_type, i, ".png")
    save_trace(result, name)
    
    # Make histogram ready list
    y <- df$y # the original observed data
    pp_se <- transform_post(y, result)
    variables = c("R", "sigma_theta", "sigma_E")
    pl_df <- pp_se [, variables]
    pl_df$sample <- 1:nrow(pl_df)
    
    if(i ==1){
      gather <- melt(pl_df, id="sample")
    }
    else{
      temp <- melt(pl_df, id="sample")
      gather =rbind(gather, temp)
    }
    
    # # extract the posterior mean from the model result
    # Rmean = summary(result)$summary["R", "mean"]
    # sigmathetamean = summary(result)$summary["sigma_theta", "mean"]
    # # save the posterior mean to the results list
    # results_list[[i]] <- c(Rmean, sigmathetamean)
  }
  
  # combine the results into a data frame
  # results_df <- do.call(rbind, results_list)
  
  # return the results
  return(gather)
}

#Run the fit several times, with same input data, but different seed for each computer model run
run_multiple_runtimes <- function(df4, df8, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, 
                              prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type, num_runs = 10) {
  #Insert with different runtimes
  # create an empty list to store the results
  results_list <- list()
  
  
  if (model_type == "w_disc") {
    stanfile <- stanfile_w_disc
  } else if (model_type == "w_E_disc") {
    stanfile <- stanfile_w_E_disc
  } else if (model_type == "wo_disc"){
    stanfile <- stanfile_wo_disc
  } else if (model_type == "p_wo_disc"){
    stanfile = stanfile_P_wo_disc
  }
  
  print(numiters)
  # loop over the desired number of runs
  for (i in 1:2) {
    # set a new random seed for each run
    seed <- 100
    if(i==1){
      df = df4
    }
    else{
      df = df8
    }
    print("HEYY")
    print(df)
    
    df2 <- create_data(df$Angle, df$Time, nc = nc, n_theta = as.integer(length(df$Time)),
                       n_E = as.integer(length(df$Time)), sigma_theta, sigma_E, prior_sigma_theta_mu,
                       prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var,
                       prior_l_delta_mu, prior_l_delta_var, prior_R_mu, prior_R_var, seed) 
    
    print("sup")
    # Create plot of df2$data_noisy_pred$t_theta and df2$data_noisy_pred$ytheta using ggplot
    p_theta <-   ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_theta, y = ytheta)) +
      geom_point() +
      xlab("Time") +
      ylab("Theta")
    plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$ytheta)
    print("yo")
    # Save plot as pdf file
    ggsave("noisy_prediction_theta.pdf", p_theta)
    
    # Create plot of df2$data_noisy_pred$t_E and df2$data_noisy_pred$yE using ggplot
    p_E <-  ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_E, y = yE)) +
      geom_point() +
      xlab("Time") +
      ylab("E")
    
    # Save plot as pdf file
    ggsave("noisy_prediction_E.pdf", p_E)
    # run the model
    result <- stan(file=stanfile,
                   data=df2$data_noisy_pred,
                   chains=3,
                   iter=numiters,
                   seed=seed)
    print(result)
    
    name = paste(model_type, i, ".png")
    save_trace(result, name)
    
    # Make histogram ready list
    y <- df$y # the original observed data
    pp_se <- transform_post(y, result)
    variables = c("R", "sigma_theta", "sigma_E")
    pl_df <- pp_se [, variables]
    pl_df$sample <- 1:nrow(pl_df)
    
    if(i ==1){
      gather <- melt(pl_df, id="sample")
    }
    else{
      temp <- melt(pl_df, id="sample")
      gather =rbind(gather, temp)
    }

  }
  return(gather)
}

#Run the fit several times, with same input data, but different seed for each computer model run
run_multiple_priors <- function(df, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, 
                                  prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type, num_runs = 10) {
  #Insert with different runtimes
  # create an empty list to store the results
  results_list <- list()
  
  
  if (model_type == "w_disc") {
    stanfile <- stanfile_w_disc
  } else if (model_type == "w_E_disc") {
    stanfile <- stanfile_w_E_disc
  } else if (model_type == "wo_disc"){
    stanfile <- stanfile_wo_disc
  } else if (model_type == "p_wo_disc"){
    stanfile = stanfile_P_wo_disc
  }
  
  print(numiters)
  # loop over the desired number of runs
  for (i in 1:3) {
    # set a new random seed for each run
    seed <- i*100
    if(i==1){
      #Informative
      prior_sigma_theta_mu = 0.03
      prior_sigma_theta_var = 0.01
      prior_R_mu = 0.2
      prior_R_var = 0.05
    }
    if(i==2){
      #Basis
      prior_sigma_theta_mu = 0
      prior_sigma_theta_var = 0.05
      prior_R_mu = 0.2
      prior_R_var = 0.1
    }
    if(i==3){
      #Weakly informative
      prior_sigma_theta_mu = 0
      prior_sigma_theta_var = 0.2
      prior_R_mu = 0.25
      prior_R_var = 0.2
    }
    print("HEYY")
    print(c(prior_sigma_theta_mu, prior_sigma_theta_var, prior_R_mu, prior_R_var))
    
    df2 <- create_data(df$Angle, df$Time, nc = nc, n_theta = as.integer(length(df$Time)),
                       n_E = as.integer(length(df$Time)), sigma_theta, sigma_E, prior_sigma_theta_mu,
                       prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var,
                       prior_l_delta_mu, prior_l_delta_var, prior_R_mu, prior_R_var, seed) 
    
    # Create plot of df2$data_noisy_pred$t_theta and df2$data_noisy_pred$ytheta using ggplot
    p_theta <-   ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_theta, y = ytheta)) +
      geom_point() +
      xlab("Time") +
      ylab("Theta")
    plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$ytheta)
    print("yo")
    # Save plot as pdf file
    ggsave("noisy_prediction_theta.pdf", p_theta)
    
    # Create plot of df2$data_noisy_pred$t_E and df2$data_noisy_pred$yE using ggplot
    p_E <-  ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_E, y = yE)) +
      geom_point() +
      xlab("Time") +
      ylab("E")
    
    # Save plot as pdf file
    ggsave("noisy_prediction_E.pdf", p_E)
    # run the model
    result <- stan(file=stanfile,
                   data=df2$data_noisy_pred,
                   chains=3,
                   iter=numiters,
                   seed=seed)
    print(result)
    
    name = paste(model_type, i, ".png")
    save_trace(result, name)
    
    # Make histogram ready list
    y <- df$y # the original observed data
    pp_se <- transform_post(y, result)
    variables = c("R", "sigma_theta", "sigma_E")
    pl_df <- pp_se [, variables]
    pl_df$sample <- 1:nrow(pl_df)
    
    if(i ==1){
      gather <- melt(pl_df, id="sample")
    }
    else{
      temp <- melt(pl_df, id="sample")
      gather =rbind(gather, temp)
    }
    
  }
  return(gather)
}

#Run the fit several times, with same input data, but different seed for each computer model run
run_different_noises_avg <- function(R_real, df, nc, theta0, L, sigma_theta_list, prior_sigma_theta_var, range_sigma_theta_upper, 
                              prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type, numruns) {
  # Define function to calculate bias and CI width
  calculate_bias_and_ci_width <- function(fit, R_real, sigma_theta_real) {
    biasR <- abs(summary(fit)$summary["R", "mean"] - R_real)
    bias_sigma_theta <- abs(summary(fit)$summary["sigma_theta", "mean"] - sigma_theta_real)
    ci_width <- summary(fit)$summary["R", "97.5%"] - summary(fit)$summary["R", "2.5%"]
    return(c(biasR, bias_sigma_theta, ci_width))
  }
  # create an empty list to store the results
  results_list <- list()
  
  
  if (model_type == "w_disc") {
    stanfile <- stanfile_w_disc
  } else if (model_type == "w_E_disc") {
    stanfile <- stanfile_w_E_disc
  } else if (model_type == "wo_disc"){
    stanfile <- stanfile_wo_disc
  }
  
  print(numiters)
  # loop over the desired number of runs
  for (i in 1:length(sigma_theta_list)) {
    # set a new random seed for each run
    bias_ci_j = 0
    for (j in 1:numruns){
      seed <- j*100
      sigma_theta = sigma_theta_list[i]
      prior_sigma_theta_mu = sigma_theta_list[i]
      df2 <- create_data(df$Angle, df$Time, nc = nc, n_theta = as.integer(length(df$Time)),
                         n_E = as.integer(length(df$Time)), sigma_theta, sigma_E, prior_sigma_theta_mu,
                         prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var,
                         prior_l_delta_mu, prior_l_delta_var, prior_R_mu, prior_R_var, seed) 
      
      print(sigma_theta)

      # Create plot of df2$data_noisy_pred$t_theta and df2$data_noisy_pred$ytheta using ggplot
      p_theta <-   ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_theta, y = ytheta)) +
        geom_point() +
        xlab("Time") +
        ylab("Theta")
      plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$ytheta)
      # Save plot as pdf file
      ggsave("noisy_prediction_theta.pdf", p_theta)
      
      # Create plot of df2$data_noisy_pred$t_E and df2$data_noisy_pred$yE using ggplot
      p_E <-  ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_E, y = yE)) +
        geom_point() +
        xlab("Time") +
        ylab("E")
      
      # Save plot as pdf file
      ggsave("noisy_prediction_E.pdf", p_E)
      # run the model
      result <- stan(file=stanfile,
                     data=df2$data_noisy_pred,
                     chains=3,
                     iter=numiters,
                     seed=seed)
      print(result)
      name = paste(model_type, i, sigma_theta, ".png")
      save_trace(result, name)
  
      
      # Make histogram ready list
      y <- df$y # the original observed data
      pp_se <- transform_post(y, result)
      variables = c("R", "sigma_theta", "sigma_E")
      pl_df <- pp_se [, variables]
      pl_df$sample <- 1:nrow(pl_df)
      
      bias_ci_j_temp = calculate_bias_and_ci_width(result, R_real, sigma_theta_list[i])
      print(bias_ci_j_temp)
      bias_ci_j = bias_ci_j + bias_ci_j_temp
      print("start if")
      print(i)
      print(j)
      }
      if(i ==1){
        print(bias_ci_j)
        bias_ci = bias_ci_j/numruns
        print("if")
        print(bias_ci)
      }
      else{
        bias_ci_temp = bias_ci_j/numruns
        bias_ci = rbind(bias_ci, bias_ci_temp)
        print("else")
        print(bias_ci)
      }
    print("finish if")
    print(bias_ci)
    
    # extract the posterior mean from the model result
    # Rmean = summary(result)$summary["R", "mean"]
    # sigmathetamean = summary(result)$summary["sigma_theta", "mean"]
    # # save the posterior mean to the results list
    # results_list[[i]] <- c(Rmean, sigmathetamean)
    
  }

  # combine the results into a data frame
  # results_df <- do.call(rbind, results_list)
  bias_ci = data.frame(R_bias = bias_ci[,1], sigma_u_bias = bias_ci[,2], CI_width = bias_ci[,3])
  # return the results
  return(bias_ci)
}




#Run the fit several times, with same input data, but different seed for each computer model run
run_different_noises <- function(R_real, df, nc, theta0, L, sigma_theta_list, prior_sigma_theta_var, range_sigma_theta_upper, 
                                 prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type, numruns) {
  # Define function to calculate bias and CI width
  calculate_bias_and_ci_width <- function(fit, R_real, sigma_theta_real) {
    biasR <- abs(summary(fit)$summary["R", "mean"] - R_real)
    bias_sigma_theta <- abs(summary(fit)$summary["sigma_theta", "mean"] - sigma_theta_real)
    ci_width <- summary(fit)$summary["R", "97.5%"] - summary(fit)$summary["R", "2.5%"]
    return(c(biasR, bias_sigma_theta, ci_width))
  }
  # create an empty list to store the results
  results_list <- list()
  
  seed = 10
  if (model_type == "w_disc") {
    stanfile <- stanfile_w_disc
  } else if (model_type == "w_E_disc") {
    stanfile <- stanfile_w_E_disc
  } else if (model_type == "wo_disc"){
    stanfile <- stanfile_wo_disc
  }
  
  print(numiters)
  # loop over the desired number of runs
  for (i in 1:length(sigma_theta_list)) {
    # set a new random seed for each run
    
    sigma_theta = sigma_theta_list[i]
    prior_sigma_theta_mu = sigma_theta_list[i]
    df2 <- create_data(df$Angle, df$Time, nc = nc, n_theta = as.integer(length(df$Time)),
                       n_E = as.integer(length(df$Time)), sigma_theta, sigma_E, prior_sigma_theta_mu,
                       prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var,
                       prior_l_delta_mu, prior_l_delta_var, prior_R_mu, prior_R_var, seed) 
    
    
      
      # Create plot of df2$data_noisy_pred$t_theta and df2$data_noisy_pred$ytheta using ggplot
      p_theta <-   ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_theta, y = ytheta)) +
        geom_point() +
        xlab("Time") +
        ylab("Theta")
      plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$ytheta)
      # Save plot as pdf file
      ggsave("noisy_prediction_theta.pdf", p_theta)
      
      # Create plot of df2$data_noisy_pred$t_E and df2$data_noisy_pred$yE using ggplot
      p_E <-  ggplot(as.data.frame(df2$data_noisy_pred), aes(x = t_E, y = yE)) +
        geom_point() +
        xlab("Time") +
        ylab("E")
      
      # Save plot as pdf file
      ggsave("noisy_prediction_E.pdf", p_E)
      # run the model
      result <- stan(file=stanfile,
                     data=df2$data_noisy_pred,
                     chains=3,
                     iter=numiters,
                     seed=seed)
      name = paste(model_type, i, sigma_theta, ".png")
      save_trace(result, name)
      
      bias_ci_temp = calculate_bias_and_ci_width(result, R_real, sigma_theta_list[i])
      # print(bias_ci_j_temp)
      # bias_ci_j = bias_ci_j + bias_ci_j_temp
      # print("start if")
      # print(i)
      # print(j)

    if(i ==1){
      bias_ci = bias_ci_temp
      print("if")
      print(bias_ci)
    }
    else{
      bias_ci = rbind(bias_ci, bias_ci_temp)
      print("else")
      print(bias_ci)
    }
  }
    print("finish if")
    print(bias_ci)
    
  bias_ci = data.frame(R_bias = bias_ci[,1], sigma_u_bias = bias_ci[,2], CI_width = bias_ci[,3])
  # return the results
  return(bias_ci)
}


mse_between_models <- function(length, process = "TP"){
  theta0list = seq(pi/10, pi/2, length.out = length)
  biaslist = rep(0, length=length(theta0list))
  for (i in 1:length(theta0list)){
    df1 = pendulum(R_real, theta0list[i], runtime, nummeasurements)
    if(process == "DP"){
      df2= pendulum_damped(L, damping, g, theta0list[i], omega0, nummeasurements, runtime, h)
    }
    else{
    df2= pendulum_true(R_real, theta0list[i], runtime, h, nummeasurements)
    }
    biaslist[i] = mean((df1$Angle-df2$Angle)**2)
  }
  df = data.frame(u_0 = theta0list, Bias = biaslist)
  return(df)
}


calculate_energy <- function(df, R, m, g, theta0) {
  L = R * g
  omega = c(0,diff(df$Angle) / diff(df$Time))
  potential = m * g * L * (1 - cos(df$Angle))
  kinetic = 0.5 * m * L^2 * omega**2
  energy = potential + kinetic
  print(potential)
  df_energy <- data.frame(Time = df$Time, Energy = energy)
  return(df_energy)
}


plot_bias_ci <- function(sigma_theta_list, result_wo_disc, result_w_E_disc, namestring) {
  
  # Plot 1: R_bias
  plot1 <- ggplot() +
    geom_line(aes(x = sigma_theta_list, y = result_wo_disc$R_bias, color = "M1_WOD"), size = 1.2, show.legend = FALSE) +
    geom_line(aes(x = sigma_theta_list, y = result_w_E_disc$R_bias, color = "M2_WDf"), size = 1.2, show.legend = FALSE) +
    xlab(expression(paste(sigma[u]))) +
    ylab("Posterior mean of R - True R") +
    ggtitle(expression(paste("Bias of R for different values of ", sigma[u])))+
    scale_color_manual(name = "Model", values = c("M1_WOD" = "#CB7370", "M2_WDf" = "#0072B2"))
  
  # Plot 2: sigma_u_bias
  plot2 <- ggplot() +
    geom_line(aes(x = sigma_theta_list, y = result_wo_disc$sigma_u_bias, color = "M1_WOD"), size = 1.2, show.legend = FALSE) +
    geom_line(aes(x = sigma_theta_list, y = results_w_E_disc$sigma_u_bias, color = "M2_WDf"), size = 1.2, show.legend = FALSE) +
    xlab(expression(paste(sigma[u]))) +
    ylab(expression(paste("Posterior mean of ", sigma[u], " - True ", sigma[u]))) +
    ggtitle(expression(paste("Bias of ", sigma[u], " for different values of ", sigma[u])))+
    scale_color_manual(name = "Model", values = c("M1_WOD" = "#CB7370", "M2_WDf" = "#0072B2"))
  
  # Plot 3: CI_width
  plot3 <- ggplot() +
    geom_line(aes(x = sigma_theta_list, y = result_wo_disc$CI_width, color = "M1_WOD"), size = 1.2) +
    geom_line(aes(x = sigma_theta_list, y = results_w_E_disc$CI_width, color = "M2_WDf"), size = 1.2) +
    xlab(expression(paste(sigma[u]))) +
    ylab("CI_width") +
    ggtitle(expression(paste("Credible Interval width of R for different ", sigma[u]))) +
    scale_color_manual(name = "Model", values = c("M1_WOD" = "#CB7370", "M2_WDf" = "#0072B2"))
  
  ggsave(paste("Rbias", namestring, ".pdf"),plot=plot1)
  ggsave(paste("sigmaubias", namestring, ".pdf"),plot=plot2)
  ggsave(paste("CIwidt", namestring, ".pdf"),plot=plot3)
  
  # Return the three plots
  return(list(plot1, plot2, plot3))
}


#Plot comparison of models with different runtimes
plot_sum_comparison_runtimes <- function(results_wo_disc, results_w_E_disc, label1, label2, R_real, sigma_theta, sigma_E,savename, numruns) {
  # extract posterior samples from both models
  num_parts <- 2
  part_size <- nrow(subset(results_wo_disc, variable == "R")) / num_parts
  # Create a blank plot
  p_wo <- ggplot() + xlab("R") + ylab("Frequency")
  p_w_E <- ggplot() + xlab("R") + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    if(i==1){
      layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD_4sec"), 
                              color = "black", alpha = 0.5)
      # Create a histogram layer for the current part with different fill color
      layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf_4sec"), 
                            color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p_wo <- p_wo + layer_wo
      p_w_E <- p_w_E + layer_w_E
    }
    if(i==2){
      layer_wo8 <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD_8sec"), 
                               color = "black", alpha = 0.5)
      # Create a histogram layer for the current part with different fill color
      layer_w_E8 <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf_8sec"), 
                                color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p_wo <- p_wo + layer_wo8
      p_w_E <- p_w_E + layer_w_E8
    }
  }
  
  p_wo  = p_wo+  geom_vline(xintercept = R_real, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD_4sec" = "lightcoral",     # reference values in data to assign colors
               "M1_WOD_8sec" = "darkred"),
    labels = c("M1_WOD_4sec" = "M1_WOD 4 sec",   
               "M1_WOD_8sec" = "M1_WOD 8 sec",         
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = "Posterior distributions of R") # Adjust the title of the fill legend
  
  p_w_E  = p_w_E+  geom_vline(xintercept = R_real, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M2_WDf_4sec" = "lightskyblue",     # reference values in data to assign colors
               "M2_WDf_8sec" = "blue4"),
    labels = c("M2_WDf_4sec" = "M2_WDf 4 sec",         
               "M2_WDf_8sec" = "M2_WDf 8 sec",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = "Posterior distributions of R") # Adjust the title of the fill legend
  
  ggsave(paste("Rposteriors_wo", savename, ".pdf"), p_wo)
  ggsave(paste("Rposteriors_w_E", savename, ".pdf"), p_w_E)

  # Create a blank plot
  p_wo <- ggplot() + xlab(expression(paste(sigma[u]))) + ylab("Frequency")
  p_w_E <- ggplot() + xlab(expression(paste(sigma[u]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "sigma_theta")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "sigma_theta")$value[(part_size * (i - 1) + 1):(part_size * i)]
    if(i==1){
      layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD_4sec"), 
                               color = "black", alpha = 0.5)
      # Create a histogram layer for the current part with different fill color
      layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf_4sec"), 
                                color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p_wo <- p_wo + layer_wo
      p_w_E <- p_w_E + layer_w_E
    }
    if(i==2){
      layer_wo8 <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD_8sec"), 
                                color = "black", alpha = 0.5)
      # Create a histogram layer for the current part with different fill color
      layer_w_E8 <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf_8sec"), 
                                 color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p_wo <- p_wo + layer_wo8
      p_w_E <- p_w_E + layer_w_E8
    }
  }
  
  p_wo  = p_wo+  geom_vline(xintercept = sigma_theta, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD_4sec" = "lightcoral",     # reference values in data to assign colors
               "M1_WOD_8sec" = "darkred"),
    labels = c("M1_WOD_4sec" = "M1_WOD 4 sec",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M1_WOD_8sec" = "M1_WOD 8 sec",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[u]))) # Adjust the title of the fill legend
  
  p_w_E  = p_w_E+  geom_vline(xintercept = sigma_theta, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M2_WDf_4sec" = "lightskyblue",     # reference values in data to assign colors
               "M2_WDf_8sec" = "blue4"),
    labels = c("M2_WDf_4sec" = "M2_WDf 4 sec",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf_8sec" = "M2_WDf 8 sec",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[u]))) # Adjust the title of the fill legend
  ggsave(paste("sigmauposteriorswo", savename, ".pdf"), p_wo)
  ggsave(paste("sigmauposteriorsw_E", savename, ".pdf"), p_w_E)
  
  
  
  # Create a blank plot
  p_wo <- ggplot() + xlab(expression(paste(sigma[f]))) + ylab("Frequency")
  p_w_E <- ggplot() + xlab(expression(paste(sigma[f]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "sigma_E")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "sigma_E")$value[(part_size * (i - 1) + 1):(part_size * i)]
    if(i==1){
      layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD_4sec"), 
                               color = "black", alpha = 0.5)
      # Create a histogram layer for the current part with different fill color
      layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf_4sec"), 
                                color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p_wo <- p_wo + layer_wo
      p_w_E <- p_w_E + layer_w_E
    }
    if(i==2){
      layer_wo8 <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD_8sec"), 
                                color = "black", alpha = 0.5)
      # Create a histogram layer for the current part with different fill color
      layer_w_E8 <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf_8sec"), 
                                 color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p_wo <- p_wo + layer_wo8
      p_w_E <- p_w_E + layer_w_E8
    }
  }
  
  p_wo  = p_wo+  geom_vline(xintercept = sigma_E, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD_4sec" = "lightcoral",     # reference values in data to assign colors
               "M1_WOD_8sec" = "darkred"),
    labels = c("M1_WOD_4sec" = "M1_WOD 4 sec",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M1_WOD_8sec" = "M1_WOD 8 sec",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[f]))) # Adjust the title of the fill legend
  
  p_w_E  = p_w_E+  geom_vline(xintercept = sigma_E, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M2_WDf_4sec" = "lightskyblue",     # reference values in data to assign colors
               "M2_WDf_8sec" = "blue4"),
    labels = c("M2_WDf_4sec" = "M2_WDf 4 sec",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf_8sec" = "M2_WDf 8 sec",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[f]))) # Adjust the title of the fill legend
  ggsave(paste("sigmafposteriorswo", savename, ".pdf"), p_wo)
  ggsave(paste("sigmafposteriorsw_E", savename, ".pdf"), p_w_E)
  
}

#Plot comparison of three models with different priors
plot_sum_comparison_diff_priors <- function(results, R_real, sigma_theta, sigma_E, savename, numruns) {
  # extract posterior samples from both models
  num_parts <- 3
  part_size <- nrow(subset(results, variable == "R")) / num_parts
  print(nrow(subset(results, variable == "R")))
  # Create a blank plot
  theme_set(
    theme_bw() + 
      theme(
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(1.8)),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 25, vjust = 0, margin= margin(b = 30))
      )
  )
  p <- ggplot() + xlab("R") + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    print(part_size)
    print(part_size * (i - 1) + 1)
    print(part_size * i)
    part_data <- subset(results, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    print(length(part_data))
    # part_data_w_E_disc <- subset(results_w_E_disc, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    if(i==1){
      layer <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Informative"), 
                               color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer
    }
    if(i==2){
      layer2 <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Basis"), 
                                color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer2
    }
    if(i==3){
      layer3 <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Weakly_informative"), 
                             color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer3
    }
  }
  
  p  = p+  geom_vline(xintercept = R_real, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("Informative" = "midnightblue",     # reference values in data to assign colors
               "Basis" = "lightblue4",
               "Weakly_informative" = "lightskyblue"),
    labels = c("Informative" = "Informative",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Basis" = "Basis",
               "Weakly_informative" = "Weakly informative",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Missing"),
    name = "Priors",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = "Posterior distributions of R") # Adjust the title of the fill legend
  ggsave(paste("Rposteriors", savename, ".pdf"), p)
  
  # Create a blank plot
  p <- ggplot() + xlab(expression(paste(sigma[u]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data<- subset(results, variable == "sigma_theta")$value[(part_size * (i - 1) + 1):(part_size * i)]
    if(i==1){
      layer <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Informative"), 
                            color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer
    }
    if(i==2){
      layer2 <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Basis"), 
                             color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer2
    }
    if(i==3){
      layer3 <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Weakly_informative"), 
                             color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer3
    }
  }
  
  p  = p+  geom_vline(xintercept = sigma_theta, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("Informative" = "midnightblue",     # reference values in data to assign colors
               "Basis" = "lightblue4",
               "Weakly_informative" = "lightskyblue"),
    labels = c("Informative" = "Informative",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Basis" = "Basis",
               "Weakly_informative" = "Weakly informative",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Missing"),
    name = "Priors",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[u]))) # Adjust the title of the fill legend
  ggsave(paste("sigmauposteriors", savename, ".pdf"), p)
  
  
  
  # Create a blank plot
  p <- ggplot() + xlab(expression(paste(sigma[f]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data <- subset(results, variable == "sigma_E")$value[(part_size * (i - 1) + 1):(part_size * i)]
    if(i==1){
      layer <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Informative"), 
                            color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer
    }
    if(i==2){
      layer2 <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Basis"), 
                             color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer2
    }
    if(i==3){
      layer3 <- geom_density(data = data.frame(x = part_data), aes(x = x, fill = "Weakly_informative"), 
                             color = "black", alpha = 0.5)
      # Add the histogram layer to the plot
      p <- p + layer3
    }
  }
  
  p  = p+  geom_vline(xintercept = sigma_E, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("Informative" = "midnightblue",     # reference values in data to assign colors
               "Basis" = "lightblue4",
               "Weakly_informative" = "lightskyblue"),
    labels = c("Informative" = "Informative",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Basis" = "Basis",
               "Weakly_informative" = "Weakly informative",          # re-label the legend (use "=" assignment to avoid mistakes)
               "Missing"),
    name = "Priors",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[f]))) # Adjust the title of the fill legend
  ggsave(paste("sigmafposteriors", savename, ".pdf"), p)
  
}

#Plot comparison of two different models
plot_sum_comparison <- function(results_wo_disc, results_w_E_disc, label1, label2, R_real, sigma_theta, sigma_E,savename, numruns) {
  # extract posterior samples from both models
  num_parts <- numruns
  part_size <- nrow(subset(results_wo_disc, variable == "R")) / num_parts
  # Create a blank plot
  p <- ggplot() + xlab("R") + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD"), 
                             color = "black", alpha = 0.5)
    # Create a histogram layer for the current part with different fill color
    layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf"), 
                              color = "black", alpha = 0.5)
    # Add the histogram layer to the plot
    p <- p + layer_wo + layer_w_E
    
  }
  
  p  = p+  geom_vline(xintercept = R_real, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD" = "lightcoral",     # reference values in data to assign colors
               "M2_WDf" = "lightskyblue"),
    labels = c("M1_WOD" = "M1_WOD",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf" = "M2_WDf",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = "Posterior distributions of R") # Adjust the title of the fill legend
  ggsave(paste("Rposteriors", savename, ".pdf"), p)
  
  
  # Create a blank plot
  p <- ggplot() + xlab(expression(paste(sigma[u]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "sigma_theta")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "sigma_theta")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD"), 
                             color = "black", alpha = 0.5)
    # Create a histogram layer for the current part with different fill color
    layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf"), 
                              color = "black", alpha = 0.5)
    # Add the histogram layer to the plot
    p <- p + layer_wo + layer_w_E
    
  }
  
  
  p  = p+  geom_vline(xintercept = sigma_theta, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD" = "lightcoral",     # reference values in data to assign colors
               "M2_WDf" = "lightskyblue"),
    labels = c("M1_WOD" = "M1_WOD",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf" = "M2_WDf",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[u]))) # Adjust the title of the fill legend
  ggsave(paste("sigmauposteriors", savename, ".pdf"), p)
  
  
  
  # Create a blank plot
  p <- ggplot() + xlab(expression(paste(sigma[f]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "sigma_E")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "sigma_E")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD"), 
                             color = "black", alpha = 0.5)
    # Create a histogram layer for the current part with different fill color
    layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf"), 
                              color = "black", alpha = 0.5)
    # Add the histogram layer to the plot
    p <- p + layer_wo + layer_w_E
    
  }
  
  
  p  = p+  geom_vline(xintercept = sigma_E, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD" = "lightcoral",     # reference values in data to assign colors
               "M2_WDf" = "lightskyblue"),
    labels = c("M1_WOD" = "M1_WOD",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf" = "M2_WDf",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[f]))) # Adjust the title of the fill legend
  ggsave(paste("sigmafposteriors", savename, ".pdf"), p)
  
}





#Plot comparison of two different models
plot_post_prior <- function(results_wo_disc, results_w_E_disc, R_real, sigma_theta, sigma_E,savename, numruns) {
  #Priors

  n <- 1000000  # Number of data points to generate
  
  # Generate normally distributed data
  Rprior <- rtruncnorm(n, a = 0, b=Inf, prior_R_mu, prior_R_var)
  sigmauprior <-rtruncnorm(n, a=0, b=Inf, prior_sigma_theta_mu, prior_sigma_theta_var)
  
  # Convert the data vector to a data frame
  df_Rprior <- data.frame(R = Rprior)
  df_sigmauprior <- data.frame(sigmau = sigmauprior)
  
  # extract posterior samples from both models
  num_parts <- numruns
  part_size <- nrow(subset(results_wo_disc, variable == "R")) / num_parts
  # Create a blank plot
  p <- ggplot() + xlab("R") + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "R")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD"), 
                             color = "black", alpha = 0.5)
    # Create a histogram layer for the current part with different fill color
    layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf"), 
                              color = "black", size=1.5, alpha = 0.5)
    # Add the histogram layer to the plot
    p <- p + layer_wo + layer_w_E
    
  }
  
  # p  = p+ layer_w_E 
  # p=p+geom_density(data=df_Rprior, aes(x=R, fill = "Prior"), alpha = 0.5) +geom_vline(xintercept = R_real, color = "springgreen4", linetype = "dashed", size=1) + layer_w_E +scale_fill_manual(                  # Manually specify fill (bar interior color)
  #   values = c("M2_WDf" = "lightskyblue",
  #              "Prior" = "orange"),
  #   labels = c("M2_WDf" = "M2_WDf",
  #              "Prior" = "Prior",
  #              "Missing"),
  #   name = "Model",                  # title of legend
  #   na.value = "grey"                 # assign a color for missing values
  # )+ ylim(0, 5)+
  #   labs(title = "Posterior distributions of R") # Adjust the title of the fill legend
  # ggsave(paste("Rposteriors", savename, ".pdf"), p)

  
  
  # 
  p  = p+ geom_density(data=df_Rprior, aes(x=R, fill = "Prior"), alpha = 0.5)
  p=p+ geom_vline(xintercept = R_real, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD" = "lightcoral",     # reference values in data to assign colors
               "M2_WDf" = "lightskyblue",
               "Prior" = "orange"),
    labels = c("M1_WOD" = "M1_WOD",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf" = "M2_WDf",
               "Prior" = "Prior",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+ 
    labs(title = "Posterior distributions of R") # Adjust the title of the fill legend
  ggsave(paste("Rposteriors", savename, ".pdf"), p)


  # Create a blank plot
  p <- ggplot() + xlab(expression(paste(sigma[u]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "sigma_theta")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "sigma_theta")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD"), 
                             color = "black", alpha = 0.5)
    # Create a histogram layer for the current part with different fill color
    layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf"), 
                              color = "black", alpha = 0.5)
    # Add the histogram layer to the plot
    p <- p + layer_wo + layer_w_E
    
  }
  
  
  p  = p+ geom_density(data=df_sigmauprior, aes(x=sigmau, fill = "Prior"), alpha = 0.5) 
  p=p+ geom_vline(xintercept = sigma_theta, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD" = "lightcoral",     # reference values in data to assign colors
               "M2_WDf" = "lightskyblue",
               "Prior" = "orange"),
    labels = c("M1_WOD" = "M1_WOD",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf" = "M2_WDf",
               "Prior" = "Prior",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[u]))) # Adjust the title of the fill legend
  ggsave(paste("sigmauposteriors", savename, ".pdf"), p)
  
  
  
  # Create a blank plot
  p <- ggplot() + xlab(expression(paste(sigma[f]))) + ylab("Frequency")
  # Iterate through each part and add the histogram layer with different colors
  for (i in 1:num_parts) {
    # Subset the data for the current part
    part_data_wo_disc <- subset(results_wo_disc, variable == "sigma_E")$value[(part_size * (i - 1) + 1):(part_size * i)]
    part_data_w_E_disc <- subset(results_w_E_disc, variable == "sigma_E")$value[(part_size * (i - 1) + 1):(part_size * i)]
    # Create a histogram layer for the current part with different fill color
    layer_wo <- geom_density(data = data.frame(x = part_data_wo_disc), aes(x = x, fill = "M1_WOD"), 
                             color = "black", alpha = 0.5)
    # Create a histogram layer for the current part with different fill color
    layer_w_E <- geom_density(data = data.frame(x = part_data_w_E_disc), aes(x = x, fill = "M2_WDf"), 
                              color = "black", alpha = 0.5)
    # Add the histogram layer to the plot
    p <- p + layer_wo + layer_w_E
    
  }
  
  
  p  = p+  geom_vline(xintercept = sigma_E, color = "springgreen4", linetype = "dashed", size=1) + scale_fill_manual(                  # Manually specify fill (bar interior color)
    values = c("M1_WOD" = "lightcoral",     # reference values in data to assign colors
               "M2_WDf" = "lightskyblue"),
    labels = c("M1_WOD" = "M1_WOD",          # re-label the legend (use "=" assignment to avoid mistakes)
               "M2_WDf" = "M2_WDf",
               "Missing"),
    name = "Model",                  # title of legend
    na.value = "grey"                 # assign a color for missing values
  )+
    labs(title = expression(paste("Posterior distributions of ", sigma[f]))) # Adjust the title of the fill legend
  ggsave(paste("sigmafposteriors", savename, ".pdf"), p)
  
}


#########################################################
#Posterior distributions and bias CI-width plots

seed = gpseed
df_LP=pendulum(R_real, theta0, runtime, nummeasurements)
df_TP=pendulum_true(R_real, theta0, runtime, h, nummeasurements)
df_DP <- pendulum_damped(L, damping, g, theta0, omega0, nummeasurements, runtime, h)
df_LPgp = simulate_pendulum_with_gp(l_delta, sigma_delta, runtime, nummeasurements, R_real, theta0, seed)


#CHOOSE SIMMODEL
if (pmodel == "LP"){
  simmodel = df_LP
  simmodel8 =pendulum(R_real, theta0, 8,2*nummeasurements)
}
if (pmodel == "LPgp"){
  simmodel = df_LPgp
  simmodel8 = simulate_pendulum_with_gp(l_delta, sigma_delta, 8, 2*nummeasurements, R_real, theta0, seed)
  df_LPgp = simulate_pendulum_with_gp(l_delta, sigma_delta, runtime, nummeasurements, R_real, theta0, seed)
}
if (pmodel == "TP"){
  simmodel = df_TP
  simmodel8 = pendulum_true(R_real,theta0, 8, h, 2*nummeasurements)
}
if (pmodel == "DP"){
  simmodel = df_DP
  simmodel8= pendulum_damped(L, damping, g, theta0, omega0, 2*nummeasurements, 8, h)
}


theme_set(
  theme_bw() + 
    theme(
      text = element_text(family = "Helvetica"),
      axis.title = element_text(size = rel(1.8)),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      plot.title = element_text(size = 25, vjust = 0, margin= margin(b = 30))
    )
)

#SUM
if(whattorun == "multipleseeds"){
  results_wo_disc = run_multiple_fits(simmodel, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "wo_disc", num_runs = numruns)
  results_w_E_disc = run_multiple_fits(simmodel, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "w_E_disc", num_runs = numruns)
  theme_set(
    theme_bw() +
      theme(
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(1.8)),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 25, vjust = 0, margin= margin(b = 30))
      )
  )
  plot_sum_comparison(results_wo_disc, results_w_E_disc, "M1_WOD", "M2_WDf", R_real, sigma_theta, sigma_E, namestring, numruns)
  }
if(whattorun == "multiplenoises"){
  #DIFFNOISES
  results_wo_disc = run_different_noises(R_real,simmodel, nc, theta0, L, sigma_theta_list, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "wo_disc", numruns = numruns)
  results_w_E_disc = run_different_noises(R_real,simmodel, nc, theta0, L, sigma_theta_list, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "w_E_disc", numruns = numruns)
  theme_set(
    theme_bw() +
      theme(
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(1.8)),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 25, vjust = 0, margin= margin(b = 30))
      )
  )
  plots = plot_bias_ci(sigma_theta_list, results_wo_disc, results_w_E_disc, namestring)

}
if(whattorun == "multipleruntimes"){
  results_wo_disc = run_multiple_runtimes(simmodel, simmodel8, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "wo_disc", num_runs = numruns)
  print(sigma_theta)
  print(prior_sigma_theta_mu)
  print(prior_sigma_theta_var)
  results_w_E_disc = run_multiple_runtimes(simmodel, simmodel8, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "w_E_disc", num_runs = numruns)
  print(namestring)
  print(results_wo_disc)
  print(results_w_E_disc)
  theme_set(
    theme_bw() +
      theme(
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(1.8)),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 25, vjust = 0, margin= margin(b = 30))
      )
  )
  plot_sum_comparison_runtimes(results_wo_disc, results_w_E_disc, "M1_WOD", "M2_WDf", R_real, sigma_theta, sigma_E, namestring, numruns)
}
if(whattorun == "multiplepriors"){
  results = run_multiple_priors(simmodel, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, 
                               prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "w_E_disc", num_runs = numruns)
  plot_sum_comparison_diff_priors(results, R_real, sigma_theta, sigma_E, namestring, numruns)
}
if(whattorun == "priorposterior"){
  results_wo_disc = run_multiple_fits(simmodel, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "wo_disc", num_runs = numruns)
  results_w_E_disc = run_multiple_fits(simmodel, nc, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, prior_sigma_delta_mu, prior_sigma_delta_var, prior_l_delta_mu, prior_l_delta_var,prior_R_mu, prior_R_var, model_type = "w_E_disc", num_runs = numruns)
  plot_post_prior(results_wo_disc, results_w_E_disc, R_real, sigma_theta, sigma_E, namestring, numruns)
}

seed = gpseed


#########################################################
#PPC

posterior_predictive <- function(pmodel, predict_stanfile_wo, predict_stanfile_w_E, nummeasurements, n_E, n_pred, seed, chains, numiters, runtime){
  if (pmodel == "LP"){
    simmodel = pendulum(R_real, theta0, runtime, nummeasurements)
  }
  if (pmodel == "LPgp"){
    simmodel = simulate_pendulum_with_gp(l_delta, sigma_delta, runtime, nummeasurements, R_real, theta0, gpseed)
  }
  if (pmodel == "TP"){
    simmodel = pendulum_true(R_real, theta0, runtime, h, nummeasurements)
  }
  if (pmodel == "DP"){
    simmodel = pendulum_damped(L, damping, g, theta0, omega0, nummeasurements, runtime, h)
  }
  if(extrapolate == "NO"){
    predictdata <- create_data_pred(simmodel$Angle,
                                    simmodel$Time,
                                    nc = 1,
                                    n_theta=as.integer(length(simmodel$Time)),
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
                                    seed = seed)
  }
  # if(extrapolate == "YES"){
  #   predictdata <- create_data_extrapolate(simmodel$Angle,
  #                                   simmodel$Time,
  #                                   nc = 1,
  #                                   n_theta=as.integer(length(simmodel$Time)),
  #                                   n_E,
  #                                   n_pred,
  #                                   sigma_theta,
  #                                   sigma_E,
  #                                   prior_sigma_theta_mu,
  #                                   prior_sigma_theta_var,
  #                                   range_sigma_theta_upper,
  #                                   prior_sigma_delta_mu,
  #                                   prior_sigma_delta_var,
  #                                   prior_l_delta_mu,
  #                                   prior_l_delta_var,
  #                                   seed = seed)
  # }

  print(pmodel)
  # plot(predictdata$data_noisy_pred$t_theta, predictdata$data_noisy_pred$ytheta)
  # points(predictdata$data_noisy_pred$t_E, predictdata$data_noisy_pred$yE, col = "blue")
  # points(predictdata$data_noisy_pred$t_theta_pred, rep(0, length=length(predictdata$data_noisy_pred$t_theta_pred)), col="RED")

  obspoints=data.frame(t_theta=predictdata$data_noisy_pred$t_theta, ytheta= predictdata$data_noisy_pred$ytheta)
  energypoints = data.frame(t_E=predictdata$data_noisy_pred$t_E, yE=predictdata$data_noisy_pred$yE)
  predpoints = data.frame(t_theta_pred=predictdata$data_noisy_pred$t_theta_pred, y_place = rep(0, length = length(predictdata$data_noisy_pred$t_theta_pred)))
  
  my_plot <- ggplot(data = obspoints, aes(x = t_theta, y = ytheta)) +
    geom_point() +
    geom_point(data = energypoints, aes(x = t_E, y = yE), col = "blue") +
    geom_point(data = predpoints, aes(x = t_theta_pred, y = y_place), col = "red")+
    ggtitle(pmodel)
  ggsave("my_plot.png", plot = my_plot, width = 6, height = 4, dpi = 300)

  
  print(predictdata)
  
  fit_Predict_wo = stan(file  = predict_stanfile_wo,
                        data = predictdata$data_noisy_pred,
                        chains =chains,
                        iter= numiters
  )
  color_scheme_set("blue")
  hist = mcmc_hist(fit_Predict_wo, pars = c("R", "sigma_theta"))
  ggsave(paste(pmodel, "histwo.png"), hist)
  acf = mcmc_acf(fit_Predict_wo, pars = c("R", "sigma_theta"))
  ggsave(paste(pmodel, "acfwo.png"), acf)
  print(fit_Predict_wo)
  color_scheme_set("mix-blue-red")
  trace = mcmc_trace(fit_Predict_wo, pars = c("R", "sigma_theta"))
  ggsave(paste(pmodel, "wo.png"), trace)
  # write.csv(fit_Predict_wo, paste(pmodel, "summarywo.txt"))
  # save_trace(fit_Predict_wo, "wo.png")
  
  fit_Predict_w_E = stan(file  = predict_stanfile_w_E,
                         data = predictdata$data_noisy_pred,
                         chains =chains,
                         iter= numiters
  )
  print(fit_Predict_w_E)
  color_scheme_set("blue")
  hist = mcmc_hist(fit_Predict_w_E, pars = c("R", "sigma_theta"))
  ggsave(paste(pmodel, "histwE.png"), hist)
  acf = mcmc_acf(fit_Predict_w_E, pars = c("R", "sigma_theta"))
  ggsave(paste(pmodel, "acfwE.png"), acf)
  color_scheme_set("mix-blue-red")
  trace = mcmc_trace(fit_Predict_w_E, pars = c("R", "sigma_theta"))
  ggsave(paste(pmodel, "wE.png"), trace)
  # write.csv(fit_Predict_w_E, paste(pmodel, "summarywo.txt"))
  # Make histogram ready list
  
  df = simmodel
  y <- df$y # the original observed data
  pp_se <- transform_post(y, fit_Predict_wo)
  variables = c("R", "sigma_theta", "sigma_E")
  pl_df <- pp_se [, variables]
  pl_df$sample <- 1:nrow(pl_df)
  
  gather_wo <- melt(pl_df, id="sample")
  
  y <- df$y # the original observed data
  pp_se <- transform_post(y, fit_Predict_w_E)
  variables = c("R", "sigma_theta", "sigma_E")
  pl_df <- pp_se [, variables]
  pl_df$sample <- 1:nrow(pl_df)
  
  gather_w_E <- melt(pl_df, id="sample")
  theme_set(
    theme_bw() +
      theme(
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(1.8)),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 25, vjust = 0, margin= margin(b = 30))
      )
  )
  
  name = paste(namestring, pmodel)
  plot_sum_comparison(gather_wo, gather_w_E, "M1_WOD", "M2_WDf", R_real, sigma_theta, sigma_E, name, numruns)
  
  return(c(predictdata,fit_wo = fit_Predict_wo, fit_w_E = fit_Predict_w_E))
}


combine_plots <- function(predictdata, pmodel, fit_Predict_wo, fit_Predict_w_E, n_pred, runtime) {
  plottitle = pmodel
  runtime_sim=runtime
  # if(extrapolate == "YES"){
  #   runtime_sim = runtime + 3
  # }
  if (pmodel == "LP"){
    simmodel = pendulum(R_real, theta0, runtime_sim, 400)
  }
  if (pmodel == "LPgp"){
    plottitle = "LP+GP"
    # if(extrapolate == "YES"){
    #   print(runtime)
    #   runtime_gp =runtime
    #   nummeasurements_gp = as.integer(400 *((runtime_gp)/runtime))
    # }
    # if(extrapolate =="NO"){
    #   nummeasurements_gp =400
    # }
    simmodel = simulate_pendulum_with_gp(l_delta, sigma_delta, runtime_sim, 400, R_real, theta0, gpseed)
    # simmodel = data.frame(head(simmodel, 400))
  }
  if (pmodel == "TP"){
    simmodel = pendulum_true(R_real, theta0, runtime_sim, h, 400)
  }
  if (pmodel == "DP"){
    simmodel = pendulum_damped(L, damping, g, theta0, omega0, 400, runtime_sim, h)
  }
  # Extract the summary table from the Stan result
  fit_summary_wo <- summary(fit_Predict_wo)
  fit_summary_w_E <- summary(fit_Predict_w_E)
  
  y_theta_summary_wo <- summary(fit_Predict_wo, pars = paste0("y_theta[", 1:n_pred, "]"))$summary
  f_theta_summary_wo <- summary(fit_Predict_wo, pars = paste0("f_theta[", 1:n_pred, "]"))$summary
  
  
  y_theta_summary_w_E <- summary(fit_Predict_w_E, pars = paste0("y_theta[", 1:n_pred, "]"))$summary
  f_theta_summary_w_E <- summary(fit_Predict_w_E, pars = paste0("f_theta[", 1:n_pred, "]"))$summary
  
  # Extract the mean values
  y_theta_means_wo <- y_theta_summary_wo[, "mean"]
  f_theta_means_wo <- f_theta_summary_wo[, "mean"]
  
  
  y_theta_means_w_E <- y_theta_summary_w_E[, "mean"]
  f_theta_means_w_E <- f_theta_summary_w_E[, "mean"]
  
  
  # Convert f_theta_means to a data frame
  y_theta_df_wo <- data.frame(t_theta = predictdata$t_theta_pred,
                              y_theta = matrix(unname(y_theta_means_wo), ncol = 1),
                              y_lower_wo = unname(y_theta_summary_wo[, "2.5%"]),
                              y_upper_wo = unname(y_theta_summary_wo[, "97.5%"]))
  # y_theta_df_wo50 <- data.frame(t_theta = predictdata$t_theta_pred,
  #                             y_theta = matrix(unname(y_theta_means_wo), ncol = 1),
  #                             y_lower_wo <- unname(y_theta_summary_wo[, "25%"]),
  #                             y_upper_wo <- unname(y_theta_summary_wo[, "75%"]))
  f_theta_df_wo <- data.frame(t_theta = predictdata$t_theta_pred,
                              f_theta = matrix(unname(f_theta_means_wo), ncol = 1))
  
  y_theta_df_w_E <- data.frame(t_theta = predictdata$t_theta_pred,
                               y_theta = matrix(unname(y_theta_means_w_E), ncol = 1),
                               y_lower_w_E = unname(y_theta_summary_w_E[, "2.5%"]),
                               y_upper_w_E = unname(y_theta_summary_w_E[, "97.5%"]))
  # y_theta_df_w_E50 <- data.frame(t_theta = predictdata$t_theta_pred,
  #                               y_theta = matrix(unname(y_theta_means_wo), ncol = 1),
  #                               y_lower_wo <- unname(y_theta_summary_wo[, "25%"]),
  #                               y_upper_wo <- unname(y_theta_summary_wo[, "75%"]))
  f_theta_df_w_E <- data.frame(t_theta = predictdata$t_theta_pred,
                               f_theta = matrix(unname(f_theta_means_w_E), ncol = 1))
  
  
  
  combined_plot <- ggplot(simmodel, aes(x = Time, y = Angle)) +
    geom_line(aes(linetype = "True"), size = 0.5) +
    geom_point(data = y_theta_df_wo, aes(x = t_theta, y = y_theta, color = "M1_WOD"), size = 0.8) +
    geom_ribbon(data = y_theta_df_wo, aes(x = t_theta, y = y_theta, ymin = y_lower_wo, ymax = y_upper_wo), fill = "#C44240", alpha = 0.2, show.legend = FALSE) +
    # geom_ribbon(data = y_theta_df_wo_50, aes(x = t_theta, y = y_theta, ymin = y_lower_wo, ymax = y_upper_wo), fill = "#C44241", alpha = 0.4, show.legend = FALSE) +
    geom_point(data = y_theta_df_w_E, aes(x = t_theta, y = y_theta, color = "M2_WDf"), size= 0.8) +
    geom_ribbon(data = y_theta_df_w_E, aes(x = t_theta, y = y_theta, ymin = y_lower_w_E, ymax = y_upper_w_E), fill = "#0072B2", alpha = 0.2, show.legend = FALSE) +
    # geom_ribbon(data = y_theta_df_w_E50, aes(x = t_theta, y = y_theta, ymin = y_lower_w_E, ymax = y_upper_w_E), fill = "#0072B3", alpha = 0.4, show.legend = FALSE) +
    geom_point(data = data.frame(t_theta = predictdata$t_theta, y_theta = predictdata$ytheta), aes(x = t_theta, y = y_theta, color = "Observed"), shape = 4, size = 2) +
    geom_line(data=simmodel, aes(x = Time, y = Angle), size = 0.5)+ggtitle(plottitle) + ylab(NULL)+xlab(NULL)+
    scale_shape_manual(values = c("True" = "black"), name = "") +
    scale_color_manual(values = c("True" = "black", "M1_WOD" = "#C44240", "M2_WDf" = "#0072B2", "Observed" = "grey0"),
                       labels = c("True" = "u(t)", "M1_WOD" = "M1_WOD", "M2_WDf" = "M2_WDf", "Observed" = "Observations"),
                       name = "") +
    scale_linetype_manual(values = c("True" = "solid"),
                          labels = c("True" = "Simulated process"),
                          name = "")
    # scale_fill_manual(values=c("95% CI M1_WOD" = "#C44240",
    #                            "95% CI M2_WDf" = "#0072B2"),
    #                   labels = c("95% CI M1_WOD"="95% CI M1_WOD",
    #                              "95% CI M2_WDf"="95% CI M2_WDf"),
    #                   name = "")
  
  if(pmodel=="LP"){
    combined_plot= combined_plot+
    labs(y = "Angle")
  }
  if(pmodel == "TP"){
    combined_plot= combined_plot+
      labs(x="Time", y = "Angle")
  }
  if(pmodel == "DP"){
    combined_plot= combined_plot+
      labs(x = "Time")
  }

  ggsave(paste(pmodel, namestring, "combined_pred.pdf"), combined_plot)
  
  return(combined_plot)
}
# # # # 

# seed =200

if(do_predictions == "YES"){
  # pmodel = "LP"
  
  theme_set(
    theme_bw() +
      theme(
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 15, vjust = 0, margin= margin(b = 10))
      )
  )

  fit_Predict_LP = posterior_predictive("LP", predict_stanfile_wo_disc, predict_stanfile_w_E_disc, nummeasurements, nummeasurements, n_pred, seed, 5, numiters, runtime)
  PPC_LP=combine_plots(fit_Predict_LP$data_noisy_pred, "LP", fit_Predict_LP$fit_wo, fit_Predict_LP$fit_w_E, n_pred, runtime)

  ggsave(paste(namestring, "LPposterior predictive.pdf"), PPC_LP, width = 10, height = 5)

  fit_Predict_LPgp = posterior_predictive("LPgp", predict_stanfile_wo_disc, predict_stanfile_w_E_disc, nummeasurements, nummeasurements, n_pred, seed, 5, numiters, runtime)
  PPC_LPgp=combine_plots(fit_Predict_LPgp$data_noisy_pred, "LPgp", fit_Predict_LPgp$fit_wo, fit_Predict_LPgp$fit_w_E, n_pred, runtime)

  ggsave(paste(namestring, "LPgpposterior predictive.pdf"), PPC_LPgp, width = 10, height = 5)

  fit_Predict_TP = posterior_predictive("TP", predict_stanfile_wo_disc, predict_stanfile_w_E_disc, nummeasurements, nummeasurements, n_pred, seed, 5, numiters, runtime)
  PPC_TP=combine_plots(fit_Predict_TP$data_noisy_pred, "TP", fit_Predict_TP$fit_wo, fit_Predict_TP$fit_w_E, n_pred, runtime)

  ggsave(paste(namestring, "TPposterior predictive.pdf"), PPC_TP, width = 10, height = 5)

  fit_Predict_DP = posterior_predictive("DP", predict_stanfile_wo_disc, predict_stanfile_w_E_disc, nummeasurements, nummeasurements, n_pred, seed, 5, numiters, runtime)
  PPC_DP=combine_plots(fit_Predict_DP$data_noisy_pred, "DP", fit_Predict_DP$fit_wo, fit_Predict_DP$fit_w_E, n_pred, runtime)

  ggsave(paste(namestring, "DPposterior predictive.pdf"), PPC_DP, width = 10, height = 5)

  # 
  # PPC_LP=combine_plots(fit_Predict_LP$data_noisy_pred, "LP", fit_Predict_LP$fit_wo, fit_Predict_LP$fit_w_E, n_pred, runtime)
  # PPC_LPgp=combine_plots(fit_Predict_LPgp$data_noisy_pred, "LPgp", fit_Predict_LPgp$fit_wo, fit_Predict_LPgp$fit_w_E, n_pred, runtime)
  # PPC_TP=combine_plots(fit_Predict_TP$data_noisy_pred, "TP", fit_Predict_TP$fit_wo, fit_Predict_TP$fit_w_E, n_pred, runtime)
  # PPC_DP=combine_plots(fit_Predict_DP$data_noisy_pred, "DP", fit_Predict_DP$fit_wo, fit_Predict_DP$fit_w_E, n_pred, runtime)

  theme_set(
    theme_bw() +
      theme(
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 15, vjust = 0, margin= margin(b = 10))
      )
  )
  
  # Plot the combined plots with the legend
  # combined_plot <- ggarrange(PPC_LP, PPC_LPgp, PPC_TP, PPC_DP, nrow = 2, ncol=2, common.legend = TRUE, legend = "bottom")
  # combined_plot
  
  # Plot the combined plots with the legend
  combined_plot <- ggarrange(
    PPC_LP, PPC_LPgp, PPC_TP, PPC_DP,
    nrow = 2, ncol = 2,
    common.legend = TRUE,
    legend = "bottom")
  
  # 
  # 
  ggsave(paste(namestring, "posterior predictive.pdf"), combined_plot, width = 10, height = 5)
}
  
