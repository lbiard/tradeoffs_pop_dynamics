# This R script can be used to simulate populations across different life histories.
# The tradeoff implemented is an intra-individual tradeoff between adult survival and fecundity

#############################################
#### Clean environment and load packages ####
#############################################

rm(list = ls())
library(ggplot2)
library(mgcv)
library(car)
library(popbio)
library(MASS)
library(psych)
library(patchwork)
library(data.table)
library(brms)
set.seed(30)

###############################
#### Create Life Histories ####
###############################

# set values for vital rates: juvenile survival, adult survival, transition rate
s1_init <- c(0.1,0.25,0.5,0.75,0.9)    # juvenile survival
s2_init <- c(0.1,0.25,0.5,0.75,0.9)    # adult survival
g_init <- c(0.1,0.25,0.5,0.75,0.9)     # transition probability from juvenile to adult

# create all combinations of vital rates
lhs <- expand.grid(s1_init, s2_init, g_init)
lhs <- as.matrix(lhs)

# solve fecundity for lambda = 1
f <- -((-1*lhs[,3]*lhs[,1]-lhs[,2]*lhs[,1]+1*lhs[,1]+lhs[,3]*lhs[,1]*lhs[,2]+1*lhs[,2]-1)/(1*lhs[,3]*lhs[,1]))    # fecundity, ensuring the asymptomatic growth rate is equal to 1

# all life histories and their associated vital rates
lhs <- cbind(lhs, f)

# Make Matrix Population Models #
vr <- cbind(lhs[,1]*(1-lhs[,3]) + lhs[,1]*lhs[,3]*lhs[,4], lhs[,1]*lhs[,3], lhs[,4]*lhs[,2], lhs[,2])

# MPM for each LH ("transition" is a list containing the transition matrices)
transition <- list()
for (i in 1:dim(vr)[1]) {
  transition[[i]] <- matrix(c(vr[i,1:4]), nrow = 2, ncol = 2)  
}

sapply(transition, lambda)  # check pop growth rate. should be the same for all MPM

# stable age distribution to use in the IBM
age_dist <- t(sapply(transition, function(x) eigen(x)$vectors[,1] / sum(eigen(x)$vectors[,1])))


# elasticity of lower level demographic parameters (using chain rule differentiation, see Caswell 2001)
sens <- t(sapply(transition, sensitivity))
s2_sens <- lhs[,4]*sens[,3] + 1*sens[,4]
s1_sens <- lhs[,3]*sens[,2] + ((lhs[,4]-1)*lhs[,3]+1)*sens[,1]
f_sens <- lhs[,2]*sens[,3] + lhs[,3]*lhs[,1]*sens[,1]
g_sens <- lhs[,1]*sens[,2] + ((lhs[,4]-1)*lhs[,1])*sens[,1]


s1_ela <- (lhs[,1]/1)*s1_sens    # elasticity of growth rate to juvenile survival
s2_ela <- (lhs[,2]/1)*s2_sens    # elasticity of growth rate to adult survival
g_ela <- (lhs[,3]/1)*g_sens      # elasticity of growth rate to maturation rate
f_ela <- (lhs[,4]/1)*f_sens      # elasticity of growth rate to fecundity

# Lower level elasticities could also be calculated with the function popbio::vitalsens, yielding same exact values

# store lower level elasticities
lower_level_elasticities <- (cbind(s1_ela, s2_ela, g_ela, f_ela))


###################################
#### Calculate Generation time ####
###################################

#### Calculated following equation (13) of Bienvenu et al. 2015, in The American Naturalist
# Since we have a mixed transition (a MPM entry that is due to both survival and fecundity),
# we have to split our transition matrice into F and S, respectively the fertility and survival matrices


# adult at time t+1 coming from an adult at time t
p22 <- lhs[,2]
# adult at time t+1 coming from juvenile at time t
p21 <- lhs[,1]*lhs[,3]
# juvenile at time t+1 coming from adult at time t
p12 <- lhs[,4]*lhs[,2]
# juvenile at time t+1 coming from juvenile at time t (SURVIVAL)
p11s <- lhs[,1]*(1-lhs[,3])
# juvenile at time t+1 coming from juvenile at time t (FECUNDITY)
p11f <- lhs[,1]*lhs[,3]*lhs[,4]

# make fertility matrix
vrf <- cbind(p11f, 0, p12, 0)
pf <- list()
for (i in 1:dim(vrf)[1]) {
  pf[[i]] <- matrix(c(vrf[i,1:4]), nrow = 2, ncol = 2)  
}

#make survival matrix
vrs <- cbind(p11s, p21, 0, p22)
ps <- list()  
for (i in 1:dim(vrs)[1]) {
  ps[[i]] <- matrix(c(vrs[i,1:4]), nrow = 2, ncol = 2)  
}

# "transition" is the full matrix (sum of pf and pm)

# to store generation time
T_gen <- list()

# equation (13) of Bienvenu et al. 2015
for (i in 1:dim(vrs)[1]) {
  T_gen[[i]] <- 1 + ((eigen(t(transition[[i]]))$vectors[,1] %*% ps[[i]]  %*% eigen(transition[[i]])$vectors[,1]) /
                       (eigen(t(transition[[i]]))$vectors[,1] %*% pf[[i]]  %*% eigen(transition[[i]])$vectors[,1]))
} 

#final baseline generation time value
gentime <- unlist(T_gen)


######################################
######################################
#### Create functions for the IBM ####
######################################
######################################


##################################################################################
#### Function to change correlation between 2 traits depending on environment ####
##################################################################################

# function to change the correlation based on the value of environment, depending on the scenario
context_dependent_cor <- function(inds, t, climate, scenario){
  
  env <- exp(climate[t])/(1+exp(climate[t])) # climate is the normally distributed environment (=logit_env), but env is bounded in 0 and 1 (see methods section)
  
  X <- cbind(inds$quality_sa, inds$quality_f)   # the 2 traits involved in the tradeoff
  
  initial_sd_1 <- sqrt(var(inds$quality_sa))  # initial sd in adult survival
  initial_sd_2 <- sqrt(var(inds$quality_f))   # initial sd in fecundity
  
  # get the eigen values/vectors
  eigen_system <- eigen(cor_matrix_tradeoff)
  
  # helper function - build rescaling matrix from eigenvalues and environment (`a`)
  mk_rescale_matrix <- function(eigen_vals, env, scenario) {
    x <- NULL
    y <- NULL
    x[1] <- 0
    y[1] <- 0
    x[2] <- sqrt(eigen_vals[1]) - 1
    y[2] <- (1 - sqrt(eigen_vals[2])) / sqrt(eigen_vals[2])
    diag(c(1 / (1 + x[scenario] * env), 1 + y[scenario] * env))
  }
  
  # 1. matrix formed of eigenvectors
  eigen_vecs <- eigen_system$vectors
  if (eigen_vecs[1,2] < 0) eigen_vecs <- -eigen_vecs
  # 2. diagonal scale matrix formed from eigenvalues and environment
  rescale_mat <- mk_rescale_matrix(eigen_system$values, env, scenario) 
  # 3. net linear transformation - depends on three matrices...
  Q <- eigen_vecs %*% rescale_mat %*% t(eigen_vecs)
  # 4. compute the context dependent values
  X2 <- X %*% Q
  # 5. rescale to original variance in each trait
  if (sqrt(var(X2)[1,1]) !=0 ){
  X2[,1] <- X2[,1]/sqrt(var(X2)[1,1]) * initial_sd_1
  }
  if (sqrt(var(X2)[2,2]) !=0 ){
  X2[,2] <- X2[,2]/sqrt(var(X2)[2,2]) * initial_sd_2
  }
   
  return(X2)       # vital rates involved in the tradeoff, with the correlation that was modified given the environment
}


############################################################
#### Function to simulate vital rates at each time step ####
############################################################

# function to get the vital rate of each individual for the given time step
simulate_mean_vr <- function(inds, t, climate, X2, lh, lower_level_elasticities){
  
  dens <- dim(inds)[1] - n_init_ind  # density negative if below initial density, so that the effect of density is positive when pop size is small, but negative when pop size is large
  logit_env <- climate[t]   # environmental value at time t
  
  # rescale amount of among-individual variance in each VR based on lower level elasticity
  gamma_survival_ju <- inds$quality_sj * sd_sj  
  gamma_transition <- inds$quality_ta * sd_ta
  gamma_survival_ad <- X2[,1] * sd_sa        # rate involved in the tradeoff
  gamma_fecundity_ad <- X2[,2] * sd_fa       # rate involved in the tradeoff
  
  ## Vital rate GLMs ##
  
  # intercepts are taken from Matrix Population Models created at the beginning
  
  # parameters for the survival juvenile glm   
  intercept_survival_ju <- logit(lh[1])
  slope_dens_survival_ju <- -0.002
  slope_env_survival_ju <- 1-lower_level_elasticities[1]    # buffering hypothesis 
  
  # parameters for the transition to adult glm   
  intercept_transition <- logit(lh[3])        
  slope_dens_transition <- -0.002
  slope_env_transition <- 1-lower_level_elasticities[3]     # buffering hypothesis 
  
  
  # parameters for the survival adult glm   
  intercept_survival_ad <- logit(lh[2])        
  slope_dens_survival_ad <- -0.002
  slope_env_survival_ad <- 1-lower_level_elasticities[2]    # buffering hypothesis 
  
  # parameters for the fecundity adult glm
  intercept_fecundity_ad <- log(lh[4])       
  slope_dens_fecundity_ad <- -0.002
  slope_env_fecundity_ad <- 1-lower_level_elasticities[4]   # buffering hypothesis 
  
  
  ## Survival juvenile GLM ##
  survival_ju <- function(env, dens) {
    
    exp(intercept_survival_ju + slope_env_survival_ju*logit_env + slope_dens_survival_ju*dens + gamma_survival_ju)/
      (1+exp(intercept_survival_ju + slope_env_survival_ju*logit_env + slope_dens_survival_ju*dens + gamma_survival_ju))
    
  }
  
  mu_sj <- survival_ju(env = logit_env, dens = dens)
  # mu_sj contain juvenile survival rate for each individual at the given time step
  
  ## Maturation GLM ##
  transition_ad <- function(env, dens) {
    
    exp(intercept_transition + slope_env_transition*logit_env + slope_dens_transition*dens + gamma_transition)/
      (1+exp(intercept_transition + slope_env_transition*logit_env + slope_dens_transition*dens + gamma_transition))
    
  }
  
  mu_ta <- transition_ad(env = logit_env, dens = dens)
  # mu_ta contain maturation rate for each individual at the given time step
  
  ## Survival adult GLM ##
  survival_ad <- function(env, dens) {
    
    exp(intercept_survival_ad + slope_env_survival_ad*logit_env + slope_dens_survival_ad*dens + gamma_survival_ad)/
      (1+exp(intercept_survival_ad + slope_env_survival_ad*logit_env + slope_dens_survival_ad*dens + gamma_survival_ad))
    
  }
  
  mu_sa <- survival_ad(env = logit_env, dens = dens)
  # mu_sa contain adult survival rate for each individual at the given time step
  
  ## Fecundity adult GLM ##
  fecundity_ad <- function(env, dens) {
    
    exp(intercept_fecundity_ad + slope_env_fecundity_ad*logit_env + slope_dens_fecundity_ad*dens + gamma_fecundity_ad)
    
  }
  
  mu_fa <- fecundity_ad(env = logit_env, dens = dens)
  # mu_fa contain fecundity rate for each individual at the given time step
  
  # store the vital rates
  vr.mu=cbind(mu_sj, mu_ta, mu_sa, mu_fa)              # vital rates per individual
  return(vr.mu)
}




##########################################################
#### Function to simulate survival process and ageing ####
##########################################################

# function to simulate the survival process, and the transition process for surviving juveniles
death <- function(inds, phi_ju, phi_ad, transition_ad){
  
  inds$stage <- ifelse(inds$stage==1 & inds$maturation==1, 2, inds$stage)
  
  #survival process
  inds$survival <- ifelse(inds$stage>1, rbinom(dim(inds)[1], 1, phi_ad[]), rbinom(dim(inds)[1], 1, phi_ju[]))
  
  #ageing process
  inds$age <- ifelse(inds$survival==1, inds$age+1, inds$age)
  
  #maturation process
  inds$maturation <- ifelse(inds$survival==1 & inds$stage==1, rbinom(dim(inds)[1], 1, transition_ad[]), 0)
  
  return(inds)
}



###################################################
#### Function to simulate reproduction process ####
###################################################

# function to simulate the reproduction process for surviving individuals
birth <- function(inds, lambda_ad, t, correlation_matrix){
  
  # surviving adults reproduce, surviving juveniles that just transitioned to adults reproduce
  inds$productivity <- ifelse(inds$survival == 1 & inds$stage > 1, rpois(dim(inds)[1], lambda_ad), 
                              ifelse(inds$survival == 1 & inds$stage == 1 & inds$maturation == 1, rpois(dim(inds)[1], lambda_ad), 0))
  
  inds$year <- t
  
  total_off <- sum(inds$productivity) # total number of offspring during the time step
  
  # create new offspring (if only 1 offspring for the population)
  if(total_off==1){
    corr.values.new <- mvrnorm(total_off, mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
    
    # ---- We now have the total number of new offspring; now add to inds
    new_inds <- data.frame(id=seq(length.out=total_off)+max(inds$id), 
                           year=rep(t, total_off), 
                           momid=rep(inds$id,inds$productivity), 
                           stage=rep(1, total_off), 
                           age=rep(0,total_off), 
                           survival=rep(1,total_off), 
                           productivity=rep(0,total_off),
                           maturation=rep(0, total_off),
                           quality_sj= corr.values.new[1], 
                           quality_ta= corr.values.new[2],
                           quality_sa= corr.values.new[3], 
                           quality_f= corr.values.new[4])
    
    # ---- Our new offspring can now be attached in the inds data frame
    inds <- rbind(inds, new_inds)
  }
  
  # create new offspring (if >1 offspring for the population)
  if(total_off>1){
    corr.values.new <- mvrnorm(total_off, mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
    
    # ---- We now have the total number of new offspring; now add to inds
    new_inds <- data.frame(id=seq(length.out=total_off)+max(inds$id),      #ID
                           year=rep(t, total_off),                         #year
                           momid=rep(inds$id,inds$productivity),           #ID of parent
                           stage=rep(1, total_off),                        #stage=1, juvenile
                           age=rep(0,total_off),                           #age 0
                           survival=rep(1,total_off),                      
                           productivity=rep(0,total_off),                  
                           maturation=rep(0, total_off),
                           quality_sj= corr.values.new[,1],                #id quality survival juvenile
                           quality_ta= corr.values.new[,2],                #id quality maturation
                           quality_sa= corr.values.new[,3],                #id quality survival adult 
                           quality_f= corr.values.new[,4])                 #id quality fecundity
    
    # ---- Our new offspring can now be attached in the inds data frame
    inds <- rbind(inds, new_inds)
  }
  
  return(inds)
}


#####################################
#####################################
#### Individual based simulation ####
#####################################
#####################################

#### Set initial conditions

time_steps <- 100    # number of time step
n_init_ind <- 1000   # initial number of individuals in the population
replicat <- 50 # number of simulation replications per life history


# 4D array to store abundance at each time step, for each simulation replicate, for each scenario, for each life history
ind_abund <- array(data = NA, dim = c(time_steps+1, replicat, 2, dim(lhs)[1]))
# 5D array to store vital rates (gamma values) at each time step, for each simulation replicate, for each scenario, for each life history
store_x <- array(data = NA, dim = c(time_steps, replicat, 2, dim(lhs)[1], 4))
# 3D array to store realized generation time for each simulation replicate, for each scenario, for each life history
gen_time <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))
# 3D array to store age at first reproduction for each simulation replicate, for each scenario, for each life history
age_first_repro_saved <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))
# 3D array to store adult survival for each simulation replicate, for each scenario, for each life history
adult_survival_saved <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))
# 3D array to store growth rate for each simulation replicate, for each scenario, for each life history
growth_rate_saved <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))


for (lh in 1:dim(lhs)[1]) {   ## loop over life histories
  
  # rescale the amount of id heterogeneity based on elasticities
  sd_sj <- 1.5-lower_level_elasticities[lh,1]     # variance in juvenile survival
  sd_ta <- 1.5-lower_level_elasticities[lh,3]     # variance in transition to adult
  sd_sa <- 1.5-lower_level_elasticities[lh,2]     # variance in adult survival
  sd_fa <- 1.5-lower_level_elasticities[lh,4]     # variance in adult fecundity
  
  for (scenario in 1:2) {    ## loop over scenarios of context dependence of covariation
    
    for (i in 1:replicat) {    ## loop over replicates of the simulation
      
      # create the environmental variable
      set.seed(i) # so that environment is different among replicates, but similar across life histories
      climate <- rnorm(time_steps, 0, 1)
      
      ts <- 0              # starting step
      
      
      ### create initial population
      inds <- as.data.frame(array(data = 0, dim = c(n_init_ind, 12)))   # initial population at time step 0
      colnames(inds) <- c("id", "year", "momid" , "stage", "age", "survival", "productivity", "maturation",
                          "quality_sj", "quality_ta", "quality_sa", "quality_f")
      inds$id <- 1:n_init_ind
      
      rho <- -0.9    # most negative correlation is -0.9
      
      correlation_matrix <- matrix(c(1, 0, 0, 0,
                                     0, 1, 0, 0,
                                     0, 0, 1, rho,
                                     0, 0, rho,1),
                                   ncol = 4, nrow = 4)
      cor_matrix_tradeoff <- correlation_matrix[3:4,3:4]
      
      # assign values to each individual for each vital rate, determine their position in the vital rate space
      corr.values <- matrix(0, nrow = n_init_ind, ncol = 4)
      corr.values <- mvrnorm(n_init_ind, mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
      
      inds$quality_sj <- corr.values[,1]  # survival juv individual heterogeneity
      inds$quality_ta <- corr.values[,2]  # transition individual heterogeneity
      inds$quality_sa <- corr.values[,3]  # survival adult individual heterogeneity
      inds$quality_f <- corr.values[,4]  # fecundity individual heterogeneity

      # initial population start with stable stage distribution
      inds$survival <- rep(1, n_init_ind)
      inds$maturation <- rep(0, n_init_ind)
      inds$stage <- c(rep(1, floor(age_dist[lh,1]*n_init_ind)), rep(2, (n_init_ind-floor(age_dist[lh,1]*n_init_ind))))
      inds$age <- rep(0, n_init_ind)         
      
      # inds_hist will store the population history of the simulation, that can be analyzed afterward
      inds_hist <- NULL
      inds_hist[[1]]  <- inds   # To store population history
      
      #### simulation, over the number of defined timesteps
      while(ts < time_steps && dim(inds[inds$survival==1,])[1]>3){
        
        inds <- subset(inds, survival==1)     # remove individuals which died last time step
        X2                 <- context_dependent_cor(inds = inds, t = ts+1, climate = climate, scenario = scenario)   # update correlation (context-dependent tradeoff) for the time step
        vr.mu              <- simulate_mean_vr(inds = inds, t = ts+1, climate = climate, X2 = X2, lh = lhs[lh,], 
                                               lower_level_elasticities = lower_level_elasticities[lh,])          # simulate vital rates for all individuals
        inds               <- death(inds, phi_ju = vr.mu[,1], phi_ad = vr.mu[,3], transition_ad = vr.mu[,2])     # survival process with maturation
        inds               <- birth(inds, lambda_ad = vr.mu[,4], t = ts+1, correlation_matrix = correlation_matrix)   # reproduction process for surviving adults
        
        ts                 <- ts + 1                # next time step
        
        inds_hist[[ts+1]]  <- inds # Add to individual history
        
        # store average vital rate values (gamma) at each time step
        store_x[ts, i, scenario, lh, 1] <- mean(inds$quality_sj)
        store_x[ts, i, scenario, lh, 2] <- mean(inds$quality_ta)
        store_x[ts, i, scenario, lh, 3] <- mean(X2[,1])
        store_x[ts, i, scenario, lh, 4] <- mean(X2[,2])
        
      }
      
      inds_hist <- rbindlist(inds_hist)   # merge all time steps into one dataset for the given replicate
      
      # extract and store abundance at each time step
      ind_abund[, i, scenario, lh] <- c(sapply(seq(0,max(inds_hist$year),1), function(x) length(inds_hist$id[inds_hist$year==x & inds_hist$survival==1])),
                                        rep(0, time_steps-max(inds_hist$year)))
      
      # age at first reproduction
      age_first_repro <- mean(subset((inds_hist[((inds_hist$survival==1 & inds_hist$stage==2 & inds_hist$age>0 & inds_hist$productivity!=0) | 
                                                   (inds_hist$survival==1 & inds_hist$stage==1 & inds_hist$maturation==1 & inds_hist$age>0 & inds_hist$productivity!=0))])[!duplicated(
                                                     (inds_hist[((inds_hist$survival==1 & inds_hist$stage==2 & inds_hist$age>0 & inds_hist$productivity!=0) |
                                                                   (inds_hist$survival==1 & inds_hist$stage==1 & inds_hist$maturation==1 & inds_hist$age>0 & inds_hist$productivity!=0))])$id)],
                                     year>50)$age)
      
      age_first_repro_saved[i, scenario, lh] <- age_first_repro
      
      # adult survival
      adult_survival <- length(inds_hist[inds_hist$age>0 & inds_hist$stage==2 & inds_hist$survival==1 & inds_hist$year>50]$survival) / 
        length(inds_hist[inds_hist$age>0 & inds_hist$stage==2 & inds_hist$year>50]$survival)
      
      adult_survival_saved[i, scenario, lh] <- adult_survival
      
      # growth rate
      growth_rate_saved[i, scenario, lh] <- geometric.mean((ind_abund[52:(time_steps+1),i,scenario,lh]/ind_abund[51:time_steps,i,scenario,lh]))
      
      # Realized generation time (set to NA if population goes extinct)
      if(growth_rate_saved[i, scenario, lh] == 0 | is.na(growth_rate_saved[i, scenario, lh])){
        gen_time[i, scenario, lh] <- NA
      } else {
        gen_time[i, scenario, lh] <- age_first_repro + (adult_survival / (geometric.mean((ind_abund[52:(time_steps+1),i,scenario,lh]/ind_abund[51:time_steps,i,scenario,lh]))
                                                                          - adult_survival))
      }
    } 
  }
  print(lh) # to have a "progress bar" during the simulation
}


#####################################################
#####################################################
#### Calculate some metrics from the simulations ####
#####################################################
#####################################################

ind_abund <- ind_abund[51:101,,,]      # discard first 50 time steps

time_steps <- 50 #amount of time steps kept


## scenario 1

# annual growth rate
stochastic_growth_rate_serie_scenario_1 <- (ind_abund[2:(time_steps+1),,1,]/ind_abund[1:time_steps,,1,])
# stochastic growth rate
stochastic_growth_rate_scenario_1 <- apply(stochastic_growth_rate_serie_scenario_1, c(2,3), function(x) geometric.mean(x[x>0]))
# coefficient of variation in annual growth rate
coef_var_sgr_scenario_1 <- apply(stochastic_growth_rate_serie_scenario_1, c(2,3), function(x)  sd(x, na.rm = T)/mean(x, na.rm = T))
# realized generation time
gen_time_scenario_1 <- gen_time[,1,]

## scenario 2

# annual growth rate
stochastic_growth_rate_serie_scenario_2 <- (ind_abund[2:(time_steps+1),,2,]/ind_abund[1:time_steps,,2,])
# stochastic growth rate
stochastic_growth_rate_scenario_2 <- apply(stochastic_growth_rate_serie_scenario_2, c(2,3), function(x) geometric.mean(x[x>0]))
# coefficient of variation in annual growth rate
coef_var_sgr_scenario_2 <- apply(stochastic_growth_rate_serie_scenario_2, c(2,3), function(x)  sd(x, na.rm = T)/mean(x, na.rm = T))
# realized generation time
gen_time_scenario_2 <- gen_time[,2,]


#####################
#### comparisons ####
#####################

# increase in context dependence, consequence for CV growth rate
change_coef_var_1 <- ((coef_var_sgr_scenario_2 - coef_var_sgr_scenario_1) / coef_var_sgr_scenario_1)*100

# increase in context dependence, consequence for realized generation time
change_gen_time_1 <- (gen_time_scenario_2 - gen_time_scenario_1) 

################################
################################
#### Figures CV growth rate ####
################################
################################

## Increase in context dependence of the TO ##

trial_plot <- as.data.frame(cbind(apply(change_coef_var_1, 2, function(x) mean(x, na.rm=T)),
                                  apply(change_coef_var_1, 2, function(x) quantile(x, probs = c(0.025), na.rm=T)),
                                  apply(change_coef_var_1, 2, function(x) quantile(x, probs = c(0.975), na.rm=T)),
                                  gentime))

point_plot <- t(change_coef_var_1)
point_plot <- c(point_plot)
point_plot <- as.data.frame(cbind(point_plot,  rep(gentime, replicat), rep((1:125), replicat)))


# change in CV growth rate ~ log(baseline generation time) + (1|life history)
fit <- brm(point_plot ~ log(V2) + (1|V3), data = point_plot, 
           cores=2, chains=2, iter = 6000, warmup = 3000)
summary(fit)

posterior <- as.matrix(fit)
dat_plot <- as.data.frame(posterior)


x2.sim <- seq(min(log(gentime)), max(log(gentime)), by = 0.005)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$b_Intercept + dat_plot$b_logV2 * (x2.sim[i])
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper)


library(grid)
text_fast <- textGrob("Fast", gp=gpar(fontsize=13, fontface="bold"))
text_slow <- textGrob("Slow", gp=gpar(fontsize=13, fontface="bold"))

p <- ggplot(trial_plot, aes(x = (x2.sim), y = bayes.c.eff.mean))+
  geom_hline(yintercept=0, linetype=2)+
  coord_cartesian(ylim = c(-50,75), clip = "off")+
  ggtitle("Increased density")+
  xlab(element_blank())+
  ylab(element_blank())+
  theme_bw() +
  theme(legend.position = "none")+
  annotation_custom(text_fast,xmin=0.3,xmax=0.3,ymin=-65,ymax=-65) + 
  annotation_custom(text_slow,xmin=2.5,xmax=2.5,ymin=-65,ymax=-65) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
p <- p + geom_point(data = point_plot, aes(y=point_plot, x=log(V2)), position = position_jitter(w = 0.1, h = 0),alpha=0.02, size=1)
p <- p + geom_line(data=plot.dat,aes(x = (x2.sim), y = bayes.c.eff.mean),color = "black", alpha = 0.8, size = 1.2)+
  geom_ribbon(data=plot.dat,aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.2)
p <- p + geom_point(data = trial_plot, aes(x=log(gentime), y=V1), size=1.5)
p
