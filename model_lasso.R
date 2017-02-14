# rm(list=ls())
require(R2jags)
require(doBy)
library(MCMCpack) # rwish function
library(plyr)
library(data.table)

########################################################
# CPE data
cpe_dat <- fread("MN_GN_CPUE_all_lakes.csv") 
head(cpe_dat)
cpe_dat[,.N] # number of obs
# Range in years
cpe_dat[,range(Year)]
# Number of lakes
cpe_dat[,length(unique(DOW))]
# Select species of interest
cpe_dat <- cpe_dat[, c("dowlknum","Year","walleye","yellow.perch","DOW"), with=FALSE]
# log-transfrom CPE
cpe_dat[, log_yp_cpe:= log(yellow.perch+0.1)]
cpe_dat[, log_wae_cpe:= log(walleye+0.1)]

# Lake predictor data
lake_dat <- fread("lake_predictors_for_ty.csv") 
head(lake_dat)
lake_dat[,.N]

# Climate data
gdd_dat <- fread("gdd_lags.csv")
head(gdd_dat)
gdd_dat[,.N]

## Process data
# Select out lakes with CPE obs n >= X
tbl1 <- table(cpe_dat$DOW)
range(tbl1)
cpe_dat <- cpe_dat[DOW %in% names(tbl1)[tbl1 >=8],]
cpe_dat[,.N]
cpe_dat[,length(unique(DOW))]
cpe_dat[,range(Year)]

# Merge GDD data with CPE data
setkey(cpe_dat, dowlknum, Year)
setkey(gdd_dat, dowlknum, Year)
cpe_gdd <- merge(cpe_dat, gdd_dat)
head(cpe_gdd)

# Number of lakes
cpe_gdd[,length(unique(DOW))]
lake_dat[,length(dowlknum)]

# Grab lake-level predictors for lakes in cpe_gdd
lake_dat_sub <- lake_dat[DOW %in% cpe_gdd$dowlknum, ]
lake_dat_sub[,length(DOW)]

# Standardize all covariates


#################################################################
########## BUGS CODE ############################################
#################################################################
# Define the model in the BUGS language and write a text file
sink("model.txt")
cat("
    model {
    # Likelihood: 
    for (i in 1:n){ 
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i]        
    } 
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha.hat[i],tau.alpha)
    mu.alpha.hat[i] <- mu.alpha
    }
    
    # Priors
    sigma ~ dunif(0, 100)
    mu.alpha ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    
    # Priors for beta parameters  
    # Bayesian LASSO -  a Laplace (double exponential) prior 
    for(k in 1:nBetas){
        b[k] ~ ddexp(0, lambda)
    }
    
    # Hyper-prior on lambda
    lambda ~ dexp(10)
    # lambda ~ dunif(0.001,10)
    # lambda~dgamma(0.1,0.1)

    # Derived quantities
    tau <- pow(sigma,-2) 
    tau.alpha <- pow(sigma.alpha,-2) 
} # end model
    ",fill = TRUE)
sink()


head(dat)
# Number of sites
J <-length(unique(dat$DOW))


dat2$G <- as.numeric(as.factor(as.numeric(dat2$DOW)))

# Load data
data <- list(y = dat$y, group = dat2$G, n = dim(dat)[1], J = J,
             x=dat$Zdd )


# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1), sigma=runif(1), sigma.alpha=runif(1) )
}


# Parameters monitored
parameters <- c("mu.alpha","alpha","sigma", "sigma.alpha","b")


# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 3


start.time = Sys.time()         # Set timer 
# Call JAGS from R 

out <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize posteriors
print(out, dig = 3)



