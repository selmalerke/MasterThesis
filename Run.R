#This code is the one executed. In this code, you choose what you want to run from "FullCode.PT"

#Disclamer. In the code, the angle u is called theta, and sigma_E=sigma_f.

namestring = "P11_DP_rt8_40_100000iters_sun" #what you want to call your run. A folder with this name is created and copies of all files needed are saved in this folder, together with the models (rds-files) and the result plots, together with diagnostics plots

#If "YES", you obtain both posterior distributions and PPCs for all processes
do_predictions = "YES"

if(do_predictions == "YES"){
  whattorun = "predictions"
}
if(do_predictions == "NO"){
  #If "NO", only posteriors distributions are runned, experiments can be choosen as follows
  #Posterior distribution, choices
  # 1. "multipleseeds" (E1) compare models on process defined in pmodel, running numruns realizations of this experiment
  # 2. "priorposterior" (E1) compare models on process defined in pmodel, plot together with prior
  # 2. "multiplenoises" (E2a) compare models for all values of sigma_u defined in sigma_theta_list, creating plots of the bias of R, the bias of sigma_u and the CI width of R
  # 3. "multipleruntimes" (E2b) compare model for observations with two different values of Delta t
  # 5. "multiplepriors" (not an experiment shown): compare models with different priors
  whattorun = "multiplenoises" #Change this to desired experiment
}

pmodel = "DP" #The process you would like to test

#If want to do extrapolation, not in the thesis
# extrapolate ="NO"

#Set seed to be able to reproduce the results
gpseed =10

# Decide parameters
L <- 2 # Length of pendulum rod. BASIS = 2
g <- 10 # Gravitational constant. BASIS = 10
R_real <- L / g # The real value, which we want to estimate BASIS = 0.2
theta0 <- pi / 5 # Starting angle in radians BASIS = pi/5
sigma_theta <- 0.03# Angle measurement noise. BASIS = 0.03
sigma_E = 0.03 #f(t) measurement noise. BASIS = 0.03
damping <- 0.5 #damping parameter for DP #BASIS = 0.5
omega0 <- 0.1  #for damped process #BASIS = 0.1
sigma_delta <- theta0^3 / 6 #for added GP in LPGP BASIS = theta0^3 / 6 
l_delta <- 1 #for added GP in LPGP BASIS = 1
h <- 0.001 #Step length for simulations of DP and TP

numiters <- 100000 #Number of MCMC iterations
sigma_theta_list = seq(0.02, 0.15, length.out = 10) #list of noises
# seeds=as.integer(seq(10, 1800, length.out = 20)) 
nummeasurements <- 40 #N_u
numruns = 1 #Number of realizations for "multipleseeds". set to 1 and 5 in E1
nc=1 #Number of replicates in observations. Set to 1
n_pred = 300 #Number of predictions in PPC
n_E = nummeasurements #N_f

#Simulation
runtime <- 8 # Delta t. Time preriod of observations
steplength <- runtime / nummeasurements # Length between each measurement


#Priors
prior_sigma_theta_mu <- 0 #BASIS = 0
prior_sigma_theta_var <- 0.02 #BASIS =0.05
range_sigma_theta_upper <- 1 # BASIS = 1
prior_sigma_delta_mu <- 0 #BASIS =0
prior_sigma_delta_var <- 1 #BASIS =1
prior_l_delta_mu <- 1 #BASIS = 1
prior_l_delta_var <- 1 #BASIS = 1
prior_R_mu <- 0.2 #BASIS = R_real
prior_R_var <- 0.1 #BASIS = 0.1

seed=28

#list of priors, for priorposterior"
prior_w_disc_mu <- c(prior_l_delta_mu, prior_sigma_delta_mu, prior_R_mu, prior_sigma_theta_mu)
prior_w_disc_var <- c(prior_l_delta_var, prior_sigma_delta_var, prior_R_var, prior_sigma_theta_var)

prior_wo_disc_mu <- c(prior_R_mu, prior_sigma_theta_mu)
prior_wo_disc_var <- c(prior_R_var, prior_sigma_theta_var)

#save true values in list
realval <- data.frame(variable = c("sigma_theta", "R", "sigma_delta", "l_delta", "sigma_E"), real = c(sigma_theta, R_real, sigma_delta, l_delta, sigma_E))

#testing = TRUE: run on computer, testing = FALSE: run on server
testing <- FALSE
# testing = TRUE

# Load functions and models from other files. Might have to change location to where your files are
if (!testing) {
  #Make new folder based on ID and copy all files there
  file = "Create_input_data.R" #used to generate noisy observations from the given process
  stanfile_wo_disc <- "C_Pendulum_model_wo_discrepancy.stan" #Model without discrepancy for estimation
  stanfile_w_disc <- "C_Pendulum_model_w_u_discrepancy.stan" #Model with discrepancy on angle for estimation
  stanfile_w_E_disc <- "C_Pendulum_model_w_E_discrepancy.stan" #Model with discrepancy on f for estimation
  # stanfile_P_wo_disc = "C_Periodic_kernel_wo_discrepancy.stan"
  predict_stanfile_wo_disc <- "C_Prediction_model_wo_discrepancy.stan" #Model without discrepancy on angle for estimation and prediction
  predict_stanfile_w_E_disc <- "C_Prediction_model_w_E_discrepancy.stan" ##Model with discrepancy on angle for estimation and prediction
  runcode = "FullCode_PT.R" #The code where all is put together and that executes the models and makes the result plots
  currentfiles = c(file, stanfile_wo_disc, stanfile_w_disc, stanfile_w_E_disc, runcode, predict_stanfile_wo_disc, predict_stanfile_w_E_disc)
  targetDirectory <- paste(namestring)
  
  # Create target directory if it doesn't exist
  if (!dir.exists(targetDirectory)) {
    dir.create(targetDirectory, recursive = TRUE)
  }
  
  # Remove spaces around '/' character in newlocation paths
  newlocation <- file.path(targetDirectory, c(file, stanfile_wo_disc, stanfile_w_disc, stanfile_w_E_disc, runcode, predict_stanfile_wo_disc, predict_stanfile_w_E_disc))
  newlocation <- gsub("\\s*/\\s*", "/", newlocation)
  
  file.copy(from=currentfiles, to=newlocation, 
            overwrite = TRUE, recursive = FALSE, 
            copy.mode = TRUE)
  source(file)
} else {
  source("~/Create_input_data.R")
  stanfile_wo_disc <- "~/C_Pendulum_model_wo_discrepancy.stan"
  stanfile_w_disc <- "~/C_Pendulum_model_w_discrepancy.stan"
  stanfile_w_E_disc <- "~/C_Pendulum_model_w_E_discrepancy.stan"
  stanfile_P_wo_disc = "~/C_periodic_kernel_wo_disc.stan"
  
  predict_stanfile_w_E_disc <- "~/C_Prediction_model_w_E_discrepancy.stan"
  predict_stanfile_wo_disc<- "~/C_Predict_model_wo_discrepancy.stan"
}


#Updates in terminal to assure that things are going as they should
print("Ready 2 run")
currentDirectory <- getwd()
print(currentDirectory)
setwd(namestring)
currentDirectory <- getwd()
print("This is the directory youre running from:")
print(currentDirectory)
print("and this is the code youre running:")
print(file.path(getwd(), runcode))

source(file.path(getwd(), runcode))
