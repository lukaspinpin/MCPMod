library(DoseFinding)
library(ggplot2)
library(logistf)
library(dplyr)
library(broom)

### logit Functions ##################################
logit <- function(p) log(p / (1 - p))
inv_logit <- function(y) 1 / (1 + exp(-y))

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

#################### Generating Data Example Function ##########################
generate_data <-function(seed=20, doses, mods, n, randMethod="PBR", B=7, x1, timetrend=FALSE, p1=0.2, pK=0.2){
  #### Data Generation ##################################
  set.seed(seed, kind = "Mersenne-Twister", sample.kind = "Rejection", normal.kind = "Inversion")
  dose_block <- c(0,10,10,25,25,100,100)
  if(randMethod=="CR"){
    dose_vector <- sample(doses, size =sum(n), replace = T, prob=n/sum(n))
  } else if (randMethod=="RA"){
    dose_vector <- rep(dose_block,sum(n)/B)
    dose_vector <- sample(dose_vector)
  } else if (randMethod=="PBR") {
    dose_vector <- rep(dose_block,sum(n)/B)
    dose_vector <- PBR(dose_vector, B)
  }
  ## emax cal
  emax_max=(logit(pK)-logit(p1))*110/100
  ## Covariate
  prob <- inv_logit(emax(dose_vector, logit(p1), emax_max, 10) + 0.6*x1 )
  ##### Time Trend ####
  if(timetrend){
    prob <- prob + seq(-0.2,0.2,length.out=sum(n)) 
  }
  prob[prob > 1] <- 1; prob[prob < 0] <- 0   # Move time trend to logit-scale
  dat <- data.frame(y = rbinom(sum(n), 1, prob), dose = dose_vector)
  return(dat)
}


#reg_mod return mean and covariance matrix or residuals
#allows for glm and logistf
#pearson and working residuals are commented out but can be used for future research
reg_mod <- function(reg="glm", covariates = FALSE, test_type="population", dat, x1, doses){
  if(test_type=="randomisation_residuals"){
    if(reg=="glm" & !covariates){
      fit <- glm(y~1 , data = dat, family = binomial)
      ### Pearson Residuals
      #residus <- residuals(fit, type="pearson")
      ### Response Residuals
      residus <- residuals(fit, type="response")
      ### Working Residuals
      #residus <- residuals(fit, type="working")
    } else if(reg=="logistf" & !covariates){
      fit <- logistf(y ~ 1, data = dat)
      ### Pearson Residuals
      #res <- dat$y - fit$predict
      #residus <- res/(sqrt(fit$hat.diag))
      ### Response Residuals
      residus <- dat$y - fit$predict
      ### Working Residuals
      #res <- dat$y - fit$predict
      #residus <- res/(fit$hat.diag)
    } else if(reg=="glm" & covariates){
      fit <- glm(y~1 + x1 , data = dat, family = binomial)
      ### Pearson Residuals
      #residus <- residuals(fit, type="pearson")
      ### Response Residuals
      residus <- residuals(fit, type="response")
      ### Working Residuals
      #residus <- residuals(fit, type="working")
    } else if(reg=="logistf" & covariates){
      fit <- logistf(y ~ 1+ x1, data = dat)
      ### Pearson Residuals
      #res <- dat$y - fit$predict
      #residus <- res/(sqrt(fit$hat.diag))
      ### Response Residuals
      residus <- dat$y - fit$predict
      ### Working Residuals
      #res <- dat$y - fit$predict
      #residus <- res/(fit$hat.diag)
    }
    return(list(residus=residus))
  } else {
    if(reg=="glm" & !covariates){
      fit <- glm(y ~ factor(dose) + 0, data = dat, family = binomial)
      mu_hat <- coef(fit); S_hat <- vcov(fit)
    } else if(reg=="logistf"  & !covariates){
      fit <- logistf(y ~ factor(dose) + 0, data = dat)
      mu_hat <- fit$coefficients; S_hat <- fit$var
    } else if(reg=="glm" & covariates){
      fit <- glm(y ~ factor(dose) + x1 -1, data = dat, family = binomial)
      mu_hat <- coef(fit)[1:length(doses)]
      S_hat <- vcov(fit)[1:length(doses),1:length(doses)]
    } else if(reg=="logistf"  & covariates){
      fit <- logistf(y ~ factor(dose) + x1 - 1, data = dat)
      mu_hat <- fit$coefficients[1:length(doses)]
      S_hat <- fit$var[1:length(doses),1:length(doses)]
    }
  }
  
  return(list(mu_hat=mu_hat, S_hat=S_hat))
}

######### MCP Test based on Person Residuals ##################################
#MCP_RD_stan <- function(dat, doses, models, contMat, residus, n) {
#  r_mean <- tapply(residus, dat$dose, mean)
 # sZ_sq <- sum((n - 1) * tapply(residus, dat$dose, var)) / ((sum(n)) - length(doses))
  
#  SZ <- apply(contMat, 2, function(col) {
 #   sum(col * r_mean) / sqrt(sZ_sq * sum(col^2 / n))
#  })
  
#  return(SZ)
#}

######### MCP Test based on Response Residuals ##################################
MCP_RD <- function(dat, doses, models, contMat, residus, n) {
  r_mean <- tapply(residus, dat$dose, mean)
  sZ_mean <- tapply(residus, dat$dose, var)
  
  SZ <- apply(contMat, 2, function(col) {
    sum(col * r_mean) / sqrt(sum(col^2 * sZ_mean / n))
  })
  
  return(SZ)
}


### rerandomizing randomisation sequence
get_randomised_data <- function(dat, randMethod, doses, B, N, n) {
  dose_vector <- switch(randMethod,
                        "CR" = sample(doses, size = N, replace = TRUE, prob = n / N),
                        "RA" = sample(dat$dose),
                        "PBR" = PBR(dat$dose, B))
  dat$dose <- dose_vector
  return(dat)
}

MCP_Rand <- function(dat=dat, doses=doses,  n=n, test_type="population", 
                     reg="glm", randMethod="PBR", B=7, x1, covariate=FALSE, nrand=1000){
  doses <- c(0, 10, 25, 100)
  mods <- DoseFinding::Mods(emax = c(10, 50), 
                            sigEmax = rbind(c(5, 3), c(25, 3)),
                            #betaMod = c(0.2475, 2.025),   # orginal betaMod
                            betaMod = c(0.1529, 0.5809),   # new betaMod
                            doses = doses)
  reg_val <- reg_mod(dat=dat, reg=reg, x1=x1, test_type=test_type, covariates=covariate, doses)
  if(test_type=="population"){
    pvals <- min(attr(MCTtest(dose=doses, resp=reg_val$mu_hat, S=reg_val$S_hat, models = mods, alternative = "one.sided", alpha =0.1, type = "general")$tStat,"pVal"))
  } else if(test_type=="randomisation_mean"){
    results <- rep(0, nrand + 1)
    results[1] <- max(MCTtest(doses, reg_val$mu_hat, S = reg_val$S_hat, models = mods, type = "general", pVal = FALSE, critV = 0)$tStat)
    for(j in 1:nrand){
      datR <- get_randomised_data(dat[, 1:2], randMethod, doses, B, sum(n), n)
      rand_reg_val <- reg_mod(dat=datR, reg=reg, x1=x1, test_type=test_type, covariates=covariate, doses)
      results[j+1] <- max(MCTtest(doses, rand_reg_val$mu_hat, S = rand_reg_val$S_hat, models = mods, type = "general", pVal = FALSE, critV = 0)$tStat)
    }
    pvals <- sum(results[-1]- results[1]>= 0 )/nrand
  } else if(test_type=="randomisation_residuals") {
    reg_val <- reg_mod(dat=dat, reg=reg, x1=x1, covariates=covariate, doses, test_type="population")
    contMat <- optContr(models = mods, w=n)$contMat
    reg_val <- reg_mod(dat=dat, reg=reg, x1=x1, covariates=covariate, test_type=test_type)
    residus <- reg_val$residus
    results <- rep(0, nrand + 1)
    results[1] <- max(MCP_RD(dat, doses, models, contMat, residus, n))
    for(j in 1:nrand){
      datR <- get_randomised_data(dat[, 1:2], randMethod, doses, B, sum(n), n)
      results[j+1] <- max(MCP_RD(datR, doses, models, contMat, residus, n))
    }
    pvals <- sum(results[-1]- results[1]>= 0 )/nrand
  }
  return(pvals)
}

