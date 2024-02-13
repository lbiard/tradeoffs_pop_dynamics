# This R script can be used to simulate populations across different life histories.
# The tradeoff implemented is an intergenerational tradeoff between parental adult survival and offspring juvenile survival

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
set.seed(3)

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



################################################################################
#### Function to create vector with given correlation to an existing vector ####
################################################################################


#### function needed to set correlation with an existing variable 
set_cor_intergen <- function(new_inds, climate, t, scenario) {
  
  env <- exp(climate[t])/(1+exp(climate[t])) # climate is the normally distributed environment (=logit_env), but env is bounded in 0 and 1 (see methods section)
  
  correlation_value <- NULL
  
  #scenario 2
  correlation_value[2] <- (env - 1)*0.9
  
  #scenario 1
  correlation_value[1] <- -0.9
  
  # make correlation matrix
  C <- matrix(correlation_value[scenario], nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  # create new vector (offspring trait)
  X2 <- rnorm(dim(new_inds)[1], 0, 1)
  # store parental and offspring trait 
  X <- cbind(new_inds$quality_sa_parent, X2)
  
  # induce correlation (does not change X[,1])
  df <- X %*% C
  
  return(df[,2])   # offspring trait, with specific correlation to parental trait
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
  gamma_survival_ad <- inds$quality_sa * sd_sa
  gamma_fecundity_ad <- inds$quality_f * sd_fa
  
  ## Vital rate GLMs ##
  
  # intercepts are taken from Matrix Population Models created at the beginning
  
  # parameters for the survival juvenile glm   
  intercept_survival_ju <- logit(lh[1])
  slope_dens_survival_ju <- -0.0002
  slope_env_survival_ju <- 1-lower_level_elasticities[1]    # buffering hypothesis 
  
  # parameters for the transition to adult glm   
  intercept_transition <- logit(lh[3])        
  slope_dens_transition <- -0.0002
  slope_env_transition <- 1-lower_level_elasticities[3]    # buffering hypothesis 
  
  
  # parameters for the survival adult glm   
  intercept_survival_ad <- logit(lh[2])        
  slope_dens_survival_ad <- -0.0002
  slope_env_survival_ad <- 1-lower_level_elasticities[2]    # buffering hypothesis 
  
  # parameters for the fecundity adult glm
  intercept_fecundity_ad <- log(lh[4])       
  slope_dens_fecundity_ad <- -0.0002
  slope_env_fecundity_ad <- 1-lower_level_elasticities[4]    # buffering hypothesis 
  
  
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
  
  vr.mu=cbind(mu_sj, mu_ta, mu_sa, mu_fa)              # vital rates per individual
  return(vr.mu)
}




##########################################################
#### Function to simulate survival process and ageing ####
##########################################################

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
birth <- function(inds, lambda_ad, t, correlation_matrix, scenario){
  
  # surviving adults reproduce, surviving juveniles that just transitioned to adults reproduce
  inds$productivity <- ifelse(inds$survival == 1 & inds$stage > 1, rpois(dim(inds)[1], lambda_ad), 
                              ifelse(inds$survival == 1 & inds$stage == 1 & inds$maturation == 1, rpois(dim(inds)[1], lambda_ad), 0)) 
  
  
  inds$year <- t
  
  total_off <- sum(inds$productivity) # total number of offspring during the time step
  
  # create new offspring (if only 1 offspring for the population)
  if(total_off==1){
    corr.values.new <- mvrnorm(total_off, mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
    
    # We have the total number of new offspring; now add to inds
    new_inds <- data.frame(id=seq(length.out=total_off)+max(inds$id), 
                           year=rep(t, total_off), 
                           momid=rep(inds$id,inds$productivity), 
                           stage=rep(1, total_off), 
                           age=rep(0,total_off), 
                           survival=rep(1,total_off), 
                           productivity=rep(0,total_off), 
                           maturation=rep(0, total_off), 
                           quality_sj = NA, 
                           quality_ta = corr.values.new[2],
                           quality_sa = corr.values.new[3], 
                           quality_f = corr.values.new[4],
                           quality_sj_parent = rep(inds$quality_sj, inds$productivity),
                           quality_ta_parent = rep(inds$quality_ta, inds$productivity),
                           quality_sa_parent = rep(inds$quality_sa, inds$productivity),
                           quality_f_parent = rep(inds$quality_f, inds$productivity))
    
    new_inds$quality_sj <- set_cor_intergen(new_inds, climate, t, scenario) # set correlation for offspring trait involved in the tradeoff
    # new offspring can now be attached in the inds data frame
    inds <- rbind(inds, new_inds)
  }
  
  # create new offspring (if >1 offspring for the population)
  if(total_off>1){
    corr.values.new <- mvrnorm(total_off, mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
    
    # We have the total number of new offspring; now add to inds
    new_inds <- data.frame(id=seq(length.out=total_off)+max(inds$id),      #ID 
                           year=rep(t, total_off),                         #year
                           momid=rep(inds$id,inds$productivity),           #ID of parent
                           stage=rep(1, total_off),                        #stage=1, juvenile
                           age=rep(0,total_off),                           #age 0
                           survival=rep(1,total_off), 
                           productivity=rep(0,total_off), 
                           maturation=rep(0, total_off), 
                           quality_sj = NA,                          #id quality survival juvenile
                           quality_ta = corr.values.new[,2],         #id quality maturation
                           quality_sa = corr.values.new[,3],         #id quality survival adult  
                           quality_f = corr.values.new[,4],          #id quality fecundity
                           quality_sj_parent = rep(inds$quality_sj, inds$productivity),
                           quality_ta_parent = rep(inds$quality_ta, inds$productivity),
                           quality_sa_parent = rep(inds$quality_sa, inds$productivity),
                           quality_f_parent = rep(inds$quality_f, inds$productivity))
    
    new_inds$quality_sj <- set_cor_intergen(new_inds, climate, t, scenario) # set correlation for offspring trait involved in the tradeoff
    # new offspring can now be attached in the inds data frame
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
# 5D array to store vital rates at each time step, for each simulation replicate, for each scenario, for each life history
store_x <- array(data = NA, dim = c(time_steps, replicat, 2, dim(lhs)[1], 4))
# 5D array to store vital rates at each time step, for each simulation replicate, for each scenario, for each life history
store_x2 <- array(data = NA, dim = c(time_steps, replicat, 2, dim(lhs)[1], 4))
# 4D array to store correlation at each time step, for each simulation replicate, for each scenario, for each life history
store_cor_x <- array(data = NA, dim = c(time_steps, replicat, 2, dim(lhs)[1]))
# 3D array to store realized generation time for each simulation replicate, for each scenario, for each life history
gen_time <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))
# 3D array to store age at first reproduction for each simulation replicate, for each scenario, for each life history
age_first_repro_saved <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))
# 3D array to store adult survival for each simulation replicate, for each scenario, for each life history
adult_survival_saved <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))
# 3D array to store growth rate for each simulation replicate, for each scenario, for each life history
growth_rate_saved <- array(data = NA, dim = c(replicat, 2, dim(lhs)[1]))


for (lh in 1:dim(lhs)[1]) {             ## loop over life histories
  
  # rescale the amount of id heterogeneity based on elasticities
  sd_sj <- 1.5-lower_level_elasticities[lh,1]     # var in juvenile survival
  sd_ta <- 1.5-lower_level_elasticities[lh,3]     # var in transition to adult
  sd_sa <- 1.5-lower_level_elasticities[lh,2]     # var in adult survival
  sd_fa <- 1.5-lower_level_elasticities[lh,4]     # var in adult fecundity
  
  for (scenario in 1:2) {      ## loop over scenarios of context dependence of covariation
    
    for (i in 1:replicat) {    ## loop over replicate of the simulation
      
      # create the environmental variable
      set.seed(i) # so that environment is different among replicates, but similar accross life histories
      climate <- rnorm(time_steps, 0, 1)
      
      ts <- 0              # starting step
      
      
      ### create initial population
      inds <- as.data.frame(array(data = 0, dim = c(n_init_ind, 16)))   # initial population at time step 0
      colnames(inds) <- c("id", "year", "momid" , "stage", "age", "survival", "productivity", "maturation",
                          "quality_sj", "quality_ta", "quality_sa", "quality_f",
                          "quality_sj_parent", "quality_ta_parent", "quality_sa_parent", "quality_f_parent")
      inds$id <- 1:n_init_ind
      
      rho <- -0.9    # most negative correlation is -0.9
      
      correlation_matrix <- matrix(c(1, 0, 0, 0, 0, 0, rho, 0,
                                     0, 1, 0, 0, 0, 0, 0, 0,
                                     0, 0, 1, 0, 0, 0, 0, 0,
                                     0, 0, 0, 1, 0, 0, 0, 0,
                                     0, 0, 0, 0, 1, 0, 0, 0,
                                     0, 0, 0, 0, 0, 1, 0, 0,
                                     rho, 0, 0, 0, 0, 0, 1, 0,
                                     0, 0, 0, 0, 0, 0, 0, 1),
                                   ncol = 8, nrow = 8)
      
      # assign values to each individual for each vital rate, determine their position in the vital rate space
      corr.values <- matrix(0, nrow = n_init_ind, ncol = 8)
      corr.values <- mvrnorm(n_init_ind, mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
      
      inds$quality_sj <- corr.values[,1]  # survival juv individual heterogeneity
      inds$quality_ta <- corr.values[,2]  # transition individual heterogeneity
      inds$quality_sa <- corr.values[,3]  # survival adult individual heterogeneity
      inds$quality_f <- corr.values[,4]  # fecundity individual heterogeneity
      
      inds$quality_sj_parent <- corr.values[,5]  # survival juv individual heterogeneity
      inds$quality_ta_parent <- corr.values[,6]  # transition individual heterogeneity
      inds$quality_sa_parent <- corr.values[,7]  # survival adult individual heterogeneity
      inds$quality_f_parent <- corr.values[,8]  # fecundity individual heterogeneity
      
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
        
        inds <- subset(inds, survival==1)      # remove individuals which died last time step
        vr.mu              <- simulate_mean_vr(inds = inds, t = ts+1, climate = climate, X2 = X2, lh = lhs[lh,], 
                                               lower_level_elasticities = lower_level_elasticities[lh,])         # simulate vital rates for all individuals
        inds               <- death(inds, phi_ju = vr.mu[,1], phi_ad = vr.mu[,3], transition_ad = vr.mu[,2])     # survival process with aging
        inds               <- birth(inds, lambda_ad = vr.mu[,4], t = ts+1, 
                                    correlation_matrix = correlation_matrix, scenario = scenario)   # reproduction process for surviving adults
        
        ts                 <- ts + 1                # next time step
        
        inds_hist[[ts+1]]  <- inds # Add to individual history
        
        # store average vital rate values (gamma) at each time step
        store_x[ts, i, scenario, lh, 1] <- mean(inds$quality_sj)
        store_x[ts, i, scenario, lh, 2] <- mean(inds$quality_ta)
        store_x[ts, i, scenario, lh, 3] <- mean(inds$quality_sa)
        store_x[ts, i, scenario, lh, 4] <- mean(inds$quality_f)
        store_cor_x[ts, i, scenario, lh] <- cor(inds$quality_sa_parent,inds$quality_sj)
        
        store_x2[ts, i, scenario, lh, 1] <- mean(vr.mu[,1])
        store_x2[ts, i, scenario, lh, 2] <- mean(vr.mu[,2])
        store_x2[ts, i, scenario, lh, 3] <- mean(vr.mu[,3])
        store_x2[ts, i, scenario, lh, 4] <- mean(vr.mu[,4])
        
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



#preprare data for regression analysis

# response variable
coef_var_1 <- tidyr::gather(as.data.frame(coef_var_sgr_scenario_1),
                            condition, cv, V1:V125)
coef_var_2 <- tidyr::gather(as.data.frame(coef_var_sgr_scenario_2),
                            condition, cv, V1:V125)
coef_var <- data.frame(coef_var_1$condition, (coef_var_2$cv - coef_var_1$cv))
colnames(coef_var) <- c("condition", "cv")

# gentime
gentime_1 <- tidyr::gather(as.data.frame(gen_time_scenario_1),
                           condition, gt, V1:V125)
gentime_2 <- tidyr::gather(as.data.frame(gen_time_scenario_2),
                           condition, gt, V1:V125)
gen_time_final <- data.frame(gentime_1$condition, (gentime_2$gt - gentime_1$gt))
colnames(gen_time_final) <- c("condition", "gt")

# mean vr
#sj
scen_2_sj <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,1], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, sj, V1:V125)
scen_1_sj <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,1], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, sj, V1:V125)
sj_final <- data.frame(scen_1_sj$condition, (scen_2_sj$sj - scen_1_sj$sj))
colnames(sj_final) <- c("condition", "sj")

#ta
scen_2_ta <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,2], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, ta, V1:V125)
scen_1_ta <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,2], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, ta, V1:V125)
ta_final <- data.frame(scen_1_ta$condition, (scen_2_ta$ta - scen_1_ta$ta))
colnames(ta_final) <- c("condition", "ta")

#sa
scen_2_sa <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,3], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, sa, V1:V125)
scen_1_sa <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,3], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, sa, V1:V125)
sa_final <- data.frame(scen_1_sa$condition, (scen_2_sa$sa - scen_1_sa$sa))
colnames(sa_final) <- c("condition", "sa")

#f
scen_2_f <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,4], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, f, V1:V125)
scen_1_f <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,4], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, f, V1:V125)
f_final <- data.frame(scen_1_f$condition, (scen_2_f$f - scen_1_f$f))
colnames(f_final) <- c("condition", "f")

# var vr

#sj
scen_2_sjv <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,1], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, sjv, V1:V125)
scen_1_sjv <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,1], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, sjv, V1:V125)
sjv_final <- data.frame(scen_1_sjv$condition, (scen_2_sjv$sjv - scen_1_sjv$sjv))
colnames(sjv_final) <- c("condition", "sjv")

#ta
scen_2_tav <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,2], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, tav, V1:V125)
scen_1_tav <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,2], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, tav, V1:V125)
tav_final <- data.frame(scen_1_tav$condition, (scen_2_tav$tav - scen_1_tav$tav))
colnames(tav_final) <- c("condition", "tav")

#sa
scen_2_sav <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,3], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, sav, V1:V125)
scen_1_sav <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,3], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, sav, V1:V125)
sav_final <- data.frame(scen_1_sav$condition, (scen_2_sav$sav - scen_1_sav$sav))
colnames(sav_final) <- c("condition", "sav")

#f
scen_2_fv <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,4], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, fv, V1:V125)
scen_1_fv <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,4], c(2,3), function(x)  (sd(x, na.rm = T)/mean(x, na.rm=T)))),
  condition, fv, V1:V125)
fv_final <- data.frame(scen_1_fv$condition, (scen_2_fv$fv - scen_1_fv$fv))
colnames(fv_final) <- c("condition", "fv")

# mean cor
cor_mean_1 <- tidyr::gather(as.data.frame(
  apply(store_cor_x[51:100,,1,], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, cor_mean, V1:V125)
cor_mean_2 <- tidyr::gather(as.data.frame(
  apply(store_cor_x[51:100,,2,], c(2,3), function(x)  mean(x, na.rm = T))),
  condition, cor_mean, V1:V125)
cor_mean <- data.frame(cor_mean_1$condition, (cor_mean_2$cor_mean - cor_mean_1$cor_mean))
colnames(cor_mean) <- c("condition", "cor_mean")

# var cor
cor_var_1 <- tidyr::gather(as.data.frame(
  apply(store_cor_x[51:100,,1,], c(2,3), function(x)  var(x, na.rm=T))),
  condition, cor_var, V1:V125)
cor_var_2 <- tidyr::gather(as.data.frame(
  apply(store_cor_x[51:100,,2,], c(2,3), function(x)  var(x, na.rm=T))),
  condition, cor_var, V1:V125)
cor_var <- data.frame(cor_var_1$condition, (cor_var_2$cor_var - cor_var_1$cor_var))
colnames(cor_var) <- c("condition", "cor_var")


# temporal correlations between variables

cor_sj_ta_1 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,c(1,2)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_sjta, V1:V125)
cor_sj_ta_2 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,c(1,2)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_sjta, V1:V125)
cor_sjta <- data.frame(cor_sj_ta_1$condition, (cor_sj_ta_2$cor_sjta - cor_sj_ta_1$cor_sjta))
colnames(cor_sjta) <- c("condition", "cor_sjta")

cor_sj_sa_1 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,c(1,3)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_sjsa, V1:V125)
cor_sj_sa_2 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,c(1,3)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_sjsa, V1:V125)
cor_sjsa <- data.frame(cor_sj_sa_1$condition, (cor_sj_sa_2$cor_sjsa - cor_sj_sa_1$cor_sjsa))
colnames(cor_sjsa) <- c("condition", "cor_sjsa")

cor_sj_f_1 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,c(1,4)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_sjf, V1:V125)
cor_sj_f_2 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,c(1,4)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_sjf, V1:V125)
cor_sjf <- data.frame(cor_sj_f_1$condition, (cor_sj_f_2$cor_sjf - cor_sj_f_1$cor_sjf))
colnames(cor_sjf) <- c("condition", "cor_sjf")

cor_ta_sa_1 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,c(2,3)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_tasa, V1:V125)
cor_ta_sa_2 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,c(2,3)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_tasa, V1:V125)
cor_tasa <- data.frame(cor_ta_sa_1$condition, (cor_ta_sa_2$cor_tasa - cor_ta_sa_1$cor_tasa))
colnames(cor_tasa) <- c("condition", "cor_tasa")

cor_ta_f_1 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,c(2,4)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_taf, V1:V125)
cor_ta_f_2 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,c(2,4)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_taf, V1:V125)
cor_taf <- data.frame(cor_ta_f_1$condition, (cor_ta_f_2$cor_taf - cor_ta_f_1$cor_taf))
colnames(cor_taf) <- c("condition", "cor_taf")

cor_sa_f_1 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,1,,c(3,4)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_saf, V1:V125)
cor_sa_f_2 <- tidyr::gather(as.data.frame(
  apply(store_x2[51:100,,2,,c(3,4)], c(2,3), function(x)  cor(x[,1], x[,2]))),
  condition, cor_saf, V1:V125)
cor_saf <- data.frame(cor_sa_f_1$condition, (cor_sa_f_2$cor_saf - cor_sa_f_1$cor_saf))
colnames(cor_saf) <- c("condition", "cor_saf")


#final data for regression analysis of delta CV growth rate vs delta mean and CV of udnerlying demographic rates
final_dataset <- data.frame(as.factor(coef_var$condition), 
                            coef_var$cv,
                            gen_time_final$gt,
                            sj_final$sj, sjv_final$sjv,
                            ta_final$ta, tav_final$tav,
                            sa_final$sa, sav_final$sav,
                            f_final$f, fv_final$fv,
                            cor_mean$cor_mean,
                            cor_var$cor_var, 
                            cor_sjta$cor_sjta, cor_sjsa$cor_sjsa, cor_sjf$cor_sjf,
                            cor_tasa$cor_tasa, cor_taf$cor_taf, cor_saf$cor_saf)

colnames(final_dataset) <- c("baseline_gt",
                             "cv",
                             "gt",
                             "sj", "sjv",
                             "ta", "tav",
                             "sa", "sav",
                             "f", "fv",
                             "cor_mean",
                             "cor_var", 
                             "cor_sjta", "cor_sjsa", "cor_sjf",
                             "cor_tasa", "cor_taf", "cor_saf")


fit <- brm(cv ~ scale(cor_mean) + scale(cor_var) +
             scale(sj) + scale(sjv) + scale(ta) + scale(tav) +
             scale(sa) + scale(sav) + scale(f) + scale(fv) +
             scale(cor_sjta) + scale(cor_sjsa) + scale(cor_sjf) +
             scale(cor_tasa) + scale(cor_taf) + scale(cor_saf) +
             (1|baseline_gt),
           data = final_dataset,
           cores=2, chains=2, iter = 4000, warmup = 2000)
summary(fit)

bayes_R2(fit)

sjPlot::tab_model(fit, 
                  pred.labels = c("Intercept", "change mean tradeoff", "change variance tradeoff", 
                                  "change mean survival juvenile", "change CV survival juvenile",
                                  "change mean maturation", "change CV maturation",
                                  "change mean survival adult", "change CV survival adult",
                                  "change mean fecundity", "change CV fecundity",
                                  "change correlation sj-ta", "change correlation sj-sa", 
                                  "change correlation sj-f", "change correlation ta-sa",
                                  "change correlation ta-f", "change correlation sa-f"),
                  dv.labels = "Tradeoff parental adult survival / offspring juvenile survival", 
                  show.icc = F, show.obs = F, show.re.var = F, digits = 3, file = "tab_sa_sj.doc")


