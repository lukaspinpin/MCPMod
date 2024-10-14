## Generate potential outcomes based on  https://doi.org/10.1002/sim.8913
## Use scenario of strong drug effect, 50 patients in total (Scenario III) and Emax model
## Code by Lukas Pin and Bjoern Bornkamp

# Load necessary libraries
library(DoseFinding)

source("MCP_Rand.R") ## load auxiliary code

export <- list(logit=logit, inv_logit=inv_logit, PBR = PBR, generate_data=generate_data, 
               reg_mod=reg_mod, MCP_RD = MCP_RD, get_randomised_data=get_randomised_data, proportion_less_than_005=proportion_less_than_005, 
               MCP_Rand=MCP_Rand, MCTtest=MCTtest, detect_separation=detect_separation, richtmyer=richtmyer, gen_pot_outcomes=gen_pot_outcomes)

gen_pot_outcomes <- function(u){
  N <- nrow(u)
  ## generate potential outcomes
  doses <- c(0, 100, 200, 400, 1000)
  ## time-point of interest
  t_end <- 12*30.4375
  ## generate individual parameters (based on quantile function)
  VA0i <- 55*exp(qnorm(u[,1], 0, sqrt(0.07)))
  VAssi <- 30*exp(qnorm(u[,2], 0, sqrt(0.2)))
  k <- 0.005*exp(qnorm(u[,3], 0, sqrt(0.2)))
  Emaxi <- 30 + qnorm(u[,4], 0, sqrt(100))
  ED50 <- 150
  eps <- qnorm(u[,5], 0, 5.3)
  Y <- matrix(nrow = N, ncol = 5)
  for(i in 1:N){
    dr <- VA0i[i] - VAssi[i]*(1 - exp(-k[i]*t_end)) + Emaxi[i]*doses/(ED50 + doses)
    Y[i,] <- dr + eps[i]
  }
  out <- data.frame(Y, base = VA0i)
  names(out)[1:5] <- paste0("Y", doses)
  out
}

#### Permuted Block Randomization ##################################
# This method is rerandomizing a sequence that was created PBR
### x is the full randomization sequence that we which to randomize
### B is the block size
PBR <- function(x, B=7){
  num_blocks <- round(length(x)/B)
  full_blocks <- num_blocks == round(length(x)/B) 
  for(i in 1:num_blocks){
    x[((i-1)*B+1):(i*B)] <- sample(x[((i-1)*B+1):(i*B)])
  }
  if(!full_blocks){
    x[num_blocks*B+1:length(x)] <- sample(x[num_blocks*B+1:length(x)])
  }
  return(x)
}

# Define the custom function
proportion_less_than_005 <- function(x) {
  # Remove NAs from the input vector
  x <- na.omit(x)
  if (length(x) == 0) return(NA) # Handle case where no data points left after omitting NAs
  sum(x < 0.05) / length(x)
}
  
  
## generate quasi-random sequence using the Richtmyer sequence
richtmyer <- function(N, s){
  primes <- c(2,3,5,7,11,13,17,19)
  if(s > length(primes))
    stop("s too large")
  gen_vec <- sqrt(primes[1:s])
  out_mat <- matrix(nrow = N, ncol = s)
  for(i in 1:N){
    out_mat[i,] <-  (i*gen_vec)%%1
  }
  out_mat
}

# Generate quasi-random sequence
n_group <- 10
u <- richtmyer(5*n_group, 5)
plot(u[,1],u[,2])
dat <- gen_pot_outcomes(u)

# Dose levels
doses <- c(0, 100, 200, 400, 1000)

## plot individual dose-response curves
matplot(doses, t(dat[,1:5]), type="l", col = 1, lty=1, ylab="VA change")
## add true population dose-response
lines(doses, colMeans(dat[,1:5]), lwd = 3, col=2)

## generate a simulated data-set, by randomly selecting from the potential outcomes
Z <- sample(rep(1:5, each=n_group)) ## random allocation rule  
resp <- dose <- numeric(5*n_group)
for(i in 1:(5*n_group)){
  dose[i] <- doses[Z[i]]
  resp[i] <- dat[i,Z[i]]
}
obs_dat <- data.frame(dose = dose, resp = resp, base = dat$base)
tapply(obs_dat$resp, obs_dat$dose, mean)
summary(lm(resp ~base, obs_dat))

# Define models for MCTtest
mods <- Mods(linear = NULL, emax = 150, sigEmax = c(175, 2),
             linlog = NULL,
             doses = doses, addArgs = list(off = 1))

# MCT Test without covariate
mm <- MCTtest(dose, resp, data=obs_dat, models = mods)
min(attr(mm$tStat, "pVal"))
#  MCT Test with adjustment for baseline
mm2 <- MCTtest(dose, resp, data=obs_dat, models = mods, addCovars = ~base)
min(attr(mm2$tStat, "pVal"))


## Scenarios
## (i) Orignial Scenario 
## (ii) Order the potential outcome data by "base" (i.e. time trend in recruitment)
## (iii) Null hypothesis: Assume all potential outcomes are equal to Y(0)

## Methods to compare
## - Population MCPMod as above (with and without covariates)
## - Randomization MCP-Mod with and without covariates (both can be based on residual based approach)

# Main function to perform pharmacological simulations
pharmaco <- function(seed=1, rand_prod = "RA", nrand =1000){
  #Parameters
  set.seed(seed=seed)
  n_group <- 10
  u <- richtmyer(5*n_group, 5)
  dat <- gen_pot_outcomes(u)
  doses <- c(0, 100, 200, 400, 1000)
  mods <- Mods(linear = NULL, emax = 150, sigEmax = c(175, 2),
               linlog = NULL,
               doses = doses, addArgs = list(off = 1))

  # Select randomization procedure
  if(rand_prod == "RA"){
    Z <- sample(rep(1:5, each=n_group)) ## random allocation rule 
  } else if (rand_prod == "PBR") {
    Z <- PBR(x=rep(c(1,1,2,2,3,3,4,4,5,5),5), B=10)   # Permuted Block Design
  }
  
   
  pvalues_all <- rep(NA,12)
  
  # Scenario (i) Original Scenario
  resp <- dose <- numeric(5*n_group)
  for(i in 1:(5*n_group)){
    dose[i] <- doses[Z[i]]
    resp[i] <- dat[i,Z[i]]
  }
  
  obs_dat <- data.frame(dose = dose, resp = resp, base = dat$base)
  
  #Population MCPMod without covariates
  mm <- MCTtest(dose, resp, data=obs_dat, models = mods)
  pvalues_all[1] = min(attr(mm$tStat, "pVal"))
  
  #Population MCPMod with covariates
  mm2 <- MCTtest(dose, resp, data=obs_dat, models = mods, addCovars = ~base)
  pvalues_all[2] = min(attr(mm2$tStat, "pVal"))
  
  #Randomization MCP-Mod without covariates 
  results <- rep(0, nrand)
  obvs <- MCTtest(dose, resp, data=obs_dat, models = mods)$tStat
  for(j in 1:nrand){
    datR <- get_randomised_data(obs_dat[, 1:2], randMethod= rand_prod, doses, B=10,  N=50, n=rep(10,5))
    results[j] <- max(MCTtest(dose, resp, data=datR, models = mods)$tStat)
  }
  pvals <- rep(0,4)
  for(i in 1:4){
    pvals[i]<- sum(results- obvs[i]>= 0 )/nrand
  }
  pvalues_all[3] = min(pvals)
  
  #Randomization MCP-Mod with covariates 
  results <- rep(0, nrand)
  obvs <- MCTtest(dose, resp, data=obs_dat, models = mods, addCovars = ~base)$tStat
  for(j in 1:nrand){
    datR <- get_randomised_data(obs_dat[, 1:3], randMethod=rand_prod, doses, B=10, N=50, n=rep(10,5))
    results[j] <- max(MCTtest(dose, resp, data=datR, models = mods, addCovars = ~base)$tStat)
  }
  pvals <- rep(0,4)
  for(i in 1:4){
    pvals[i]<- sum(results- obvs[i]>= 0 )/nrand
  }
  pvalues_all[4] = min(pvals)
  
  # Scenario (ii) Order the potential outcome data by "base"
  dat_T <- dat[order(dat$base), ]
  resp <- dose <- numeric(5*n_group)
  for(i in 1:(5*n_group)){
    dose[i] <- doses[Z[i]]
    resp[i] <- dat_T[i,Z[i]]
  }
  obs_dat <- data.frame(dose = dose, resp = resp, base = dat_T$base)
  
  #Population MCPMod without covariates
  mm <- MCTtest(dose, resp, data=obs_dat, models = mods)
  pvalues_all[5] = min(attr(mm$tStat, "pVal"))
  
  #Population MCPMod with covariates
  mm2 <- MCTtest(dose, resp, data=obs_dat, models = mods, addCovars = ~base)
  pvalues_all[6] = min(attr(mm2$tStat, "pVal"))
  
  #Randomization MCP-Mod without covariates 
  results <- rep(0, nrand)
  obvs <- MCTtest(dose, resp, data=obs_dat, models = mods)$tStat
  for(j in 1:nrand){
    datR <- get_randomised_data(obs_dat[, 1:2], randMethod=rand_prod, doses, B=10, N=50, n=rep(10,5))
    results[j] <- max(MCTtest(dose, resp, data=datR, models = mods)$tStat)
  }
  pvals <- rep(0,4)
  for(i in 1:4){
    pvals[i]<- sum(results- obvs[i]>= 0 )/nrand
  }
  pvalues_all[7] = min(pvals)
  
  #Randomization MCP-Mod with covariates 
  results <- rep(0, nrand)
  obvs <- MCTtest(dose, resp, data=obs_dat, models = mods, addCovars = ~base)$tStat
  for(j in 1:nrand){
    datR <- get_randomised_data(obs_dat[, 1:3], randMethod=rand_prod, doses, B=10, N=50, n=rep(10,5))
    results[j] <- max(MCTtest(dose, resp, data=datR, models = mods, addCovars = ~base)$tStat)
  }
  pvals <- rep(0,4)
  for(i in 1:4){
    pvals[i]<- sum(results- obvs[i]>= 0 )/nrand
  }
  pvalues_all[8] = min(pvals)
  
  ## (iii) Null hypothesis: Assume all potential outcomes are equal to Y(0)
  dat_N <- dat
  dat_N[,2:5] <- dat_N[,1] 
  resp <- dose <- numeric(5*n_group)
  for(i in 1:(5*n_group)){
    dose[i] <- doses[Z[i]]
    resp[i] <- dat_N[i,Z[i]]
  }
  obs_dat <- data.frame(dose = dose, resp = resp, base = dat_N$base)
  
  #Population MCPMod without covariates
  mm <- MCTtest(dose, resp, data=obs_dat, models = mods)
  pvalues_all[9] = min(attr(mm$tStat, "pVal"))
  
  #Population MCPMod with covariates
  mm2 <- MCTtest(dose, resp, data=obs_dat, models = mods, addCovars = ~base)
  pvalues_all[10] = min(attr(mm2$tStat, "pVal"))
  results <- rep(0, nrand)
  
  #Randomization MCP-Mod without covariates 
  obvs <- MCTtest(dose, resp, data=obs_dat, models = mods)$tStat
  for(j in 1:nrand){
    datR <- get_randomised_data(obs_dat[, 1:2], randMethod=rand_prod, doses, B=10, N=50, n=rep(10,5))
    results[j] <- max(MCTtest(dose, resp, data=datR, models = mods)$tStat)
  }
  pvals <- rep(0,4)
  for(i in 1:4){
    pvals[i]<- sum(results- obvs[i]>= 0 )/nrand
  }
  pvalues_all[11] = min(pvals)
  
  #Randomization MCP-Mod with covariates 
  results <- rep(0, nrand)
  obvs <- MCTtest(dose, resp, data=obs_dat, models = mods, addCovars = ~base)$tStat
  for(j in 1:nrand){
    datR <- get_randomised_data(obs_dat[, 1:3], randMethod=rand_prod, doses, B=10, N=50, n=rep(10,5))
    results[j] <- max(MCTtest(dose, resp, data=datR, models = mods, addCovars = ~base)$tStat)
  }
  pvals <- rep(0,4)
  for(i in 1:4){
    pvals[i]<- sum(results- obvs[i]>= 0 )/nrand
  }
  pvalues_all[12] = min(pvals)
  return(pvalues_all)
}

grid <- expand.grid(
  rand_prod = c("RA", "PBR") ,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE)


for (i in 1:nrow(grid)) {
  param_constatnt <- list(nrand=1000, rand_prod= grid[i,])
  
  nsim = 10000
  seeds <- data.frame(seed=seq(1,nsim))

  sim_out <- Q_rows(fun = pharmaco, 
                    df = seeds,
                    const = param_constatnt,
                    seed = 23456,
                    n_jobs = nrow(seeds),
                    template = list(
                      walltime = 300,
                      job_name = "Pharmaco",
                      log_file = "logs.txt"),
                    pkgs = c("data.table", "DoseFinding", "logistf"),
                    export = export)
  
  filename <- paste("simulation_pharmaco", grid[i,], param_constatnt$nrand, "FULL", sep="_")
  saveRDS(sim_out, file= filename)
  
  sim_matrix_p_vals <- as.data.frame(t(do.call(rbind, sim_out)))
  
  results <- cbind(TypeIerror_Power=apply(sim_matrix_p_vals, 1, proportion_less_than_005))
  print(results)
  
  filename <- paste("simulation_pharmaco", grid[i,], param_constatnt$nrand, "RESULTS", sep="_")
  saveRDS(results, file= filename)
}