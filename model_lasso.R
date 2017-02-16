# rm(list=ls())
require(R2jags)
library(MCMCpack) # rwish function
library(plyr)
library(data.table)
library(car)

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
lake_dat_sub <- lake_dat[DOW %in% cpe_gdd$DOW, ]
lake_dat_sub[,length(DOW)]
lake_dat_sub[,.N]
# write.csv(lake_dat_sub,"lake_data.csv",row.names=F)
summary(lake_dat_sub)

# Standardize all covariates
cpe_gdd[, z_wae_cpe:= as.numeric(scale(log_wae_cpe))]
cpe_gdd[, z_gdd:= as.numeric(scale(gdd))]
cpe_gdd[, z_gdd1:= as.numeric(scale(gdd.lag1))]
cpe_gdd[, z_gdd2:= as.numeric(scale(gdd.lag2))]
cpe_gdd[, z_gdd3:= as.numeric(scale(gdd.lag3))]
cpe_gdd[, z_gdd4:= as.numeric(scale(gdd.lag4))]
cpe_gdd[, z_gdd5:= as.numeric(scale(gdd.lag5))]
cpe_gdd[, z_gddma:= as.numeric(scale(gdd.ma.5))]

head(cpe_gdd)
# Correlation among gdd predictors and wae
gdd_cor <- cpe_gdd[,c(15:22)]
cor(gdd_cor)

# Transform and standardize lake predictors
lake_dat_sub[, z_area := as.numeric(scale(log(area.hectares)))]
lake_dat_sub[, z_littoral := as.numeric(scale(prop.littoral))]
lake_dat_sub[, z_depth := as.numeric(scale(log(max_depth_m)))]
lake_dat_sub[, z_secchi := as.numeric(scale(Secchi.lake.mean))]
lake_dat_sub[, z_alk := as.numeric(scale(alkalinity))]
lake_dat_sub[, z_gddMean := as.numeric(scale(mean.gdd))]

lake_cor <- lake_dat_sub[,c(9:14)]
cor(lake_cor,use = "complete.obs")

# hist(lake_dat_sub$alkalinity)
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
      mu[i] <- alpha[group[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
                                b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] + b[8] * x8[i]       
    } 
    
    # Level-2 of the model
    for(j in 1:J){
      alpha[j] ~ dnorm(mu.alpha.hat[j],tau.alpha)
      mu.alpha.hat[j] <- mu.alpha + s[1] * z1[j] + s[2] * z2[j] + s[3] * z3[j] + s[4] * z4[j] 
                          # + s[5] * z5[j]+ s[6] * z6[j] 
    }
    
    # Priors
    sigma ~ dunif(0, 100)
    mu.alpha ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    
    # Priors for beta parameters  
    # Bayesian LASSO -  a Laplace (double exponential) prior
    # Level-1 predictors
    for(k in 1:8){
        b[k] ~ ddexp(0, lambda)
    }

    # Level-2 predictors
    for(k in 1:4){
        s[k] ~ ddexp(0, lambda2)
    }
    
    # Hyper-prior on lambda
    # lambda ~ dexp(10)
    lambda ~ dgamma(0.1,0.1)
    lambda2 ~ dgamma(0.1,0.1)
    # lambda~dgamma(0.01,0.01)

    # Derived quantities
    tau <- pow(sigma,-2) 
    tau.alpha <- pow(sigma.alpha,-2) 

# Calculate marginal and conditional R2
  for (i in 1:n){
    predictY[i] <- mu[i]

}

 Vfixed<-(sd(predictY))^2
 Vresidual <- 1/tau # get the variance of residuals
 Vrandom <- 1/tau.alpha  # get the variance of random intercept
 marginalR2 <- Vfixed/(Vfixed+Vrandom+Vresidual) # calculate marginalR2 (fixed effects only)
 conditionalR2 <- (Vrandom+Vfixed)/(Vfixed+Vrandom+Vresidual) # calculate conditional R2 (fixed + random)

} # end model
    ",fill = TRUE)
sink()


head(cpe_gdd)
# Number of sites
J <-length(unique(cpe_gdd$DOW))

# Site indicator
cpe_gdd[, G := as.numeric(as.factor(as.numeric(DOW)))]

# Load data
data <- list(y = cpe_gdd$log_yp_cpe, group = cpe_gdd$G, n = cpe_gdd[,.N], J = J,
             x1 = cpe_gdd$z_gdd, x2 = cpe_gdd$z_gdd1, x3 = cpe_gdd$z_gdd2, x4 = cpe_gdd$z_gdd3,
             x5 = cpe_gdd$z_gdd4, x6 = cpe_gdd$z_gdd5, x7 = cpe_gdd$z_gddma, x8 = cpe_gdd$z_wae_cpe,
             z1 = lake_dat_sub$z_area, z2 = lake_dat_sub$z_littoral, z3 = lake_dat_sub$z_depth, z4 = lake_dat_sub$z_gddMean)


# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1), sigma=runif(1), sigma.alpha=runif(1) )
}


# Parameters monitored
parameters <- c("mu.alpha","sigma", "sigma.alpha","b","lambda",
                "predictY","marginalR2","conditionalR2","lambda2","s")


# MCMC settings
ni <- 90000
nt <- 2
nb <- 70000
nc <- 3


start.time = Sys.time()         # Set timer 
# Call JAGS from R 

out <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time

# Find which parameters, if any, have Rhat > 1.1
which(out$BUGSoutput$summary[, c("Rhat")] > 1.1)


# Summarize posteriors
print(out, dig = 3)
# traceplot(out)

# str(out)

betaEsts <- matrix(NA, nrow=3,ncol=8)
for(i in 1:8){
  betaEsts[1,i] <- mean(quantile(out$BUGSoutput$sims.list$b[,i],c(0.05,0.95)))
  betaEsts[2:3,i] <- quantile(out$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
betaEsts


sEsts <- matrix(NA, nrow=3,ncol=4)
for(i in 1:4){
  sEsts[1,i] <- mean(quantile(out$BUGSoutput$sims.list$s[,i],c(0.05,0.95)))
  sEsts[2:3,i] <- quantile(out$BUGSoutput$sims.list$s[,i],c(0.05,0.95))
}
sEsts


predicted <- out$BUGSoutput$mean$predictY
observed <- cpe_gdd$log_yp_cpe

res <- 6
name_figure <- "observed_predicted_cpe.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 		# save default, for resetting...

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.5,0.5,0,0), oma=c(2,2,0,0) )

size.labels = 1
size.text = 1.0
y.label = expression(paste('Predicted ',log[e],'(CPE)' ))
x.label = expression(paste('Observed ',log[e],'(CPE)' ))

plot(predicted ~ observed, type='n', axes=F, ylim=c(min(predicted), max(predicted)),
     xlab='',ylab='', xlim=c(min(observed),max(observed)) )
axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0.3,0) )
axis(side=2,cex.axis=size.text, las=1, mgp=c(0,0.3,0),tck=-0.01)

# Add data points
points(observed, predicted, pch=16, cex=1.0, col='black')

# Add axis labels
mtext(x.label, line = 0.7, side = 1, cex = size.text, outer=T)
mtext(y.label, line = 0.5, side = 2, cex = size.text, outer=T)
abline(0,1, lwd=2, col="blue")

box()
par(def.par)
dev.off()
### END PLOT

# out.array <- out$BUGSoutput$sims.array
# output <- rbind(out.array[,1,],out.array[,2,],out.array[,3,]) #convert to matrix; the matrix directly from out output (out$BUGSoutput$sims.matrix) does not sort by chain
# 
# # Number of samples retained
# n.keep <-out$BUGSoutput$n.keep
# 
# # Trace plots
# y.min=min(output[,'lambda'],na.rm=T)
# y.max=max(output[,'lambda'],na.rm=T)
# plot(x=1:n.keep,y=output[1:n.keep,'lambda'],type='l',lwd=1,col='red',ylim=c(y.min,y.max),main=NULL)
# lines(x=1:n.keep,y=output[(n.keep+1):(n.keep*2),'lambda'],type='l',lwd=1,col='green')
# lines(x=1:n.keep,y=output[(2*n.keep+1):(n.keep*3),'lambda'],type='l',lwd=1,col='blue') #3 chains in total
# 
# y.min=min(output[,'b[1]'],na.rm=T)
# y.max=max(output[,'b[1]'],na.rm=T)
# plot(x=1:n.keep,y=output[1:n.keep,'b[1]'],type='l',lwd=1,col='red',ylim=c(y.min,y.max),main=NULL)
# lines(x=1:n.keep,y=output[(n.keep+1):(n.keep*2),'b[1]'],type='l',lwd=1,col='green')
# lines(x=1:n.keep,y=output[(2*n.keep+1):(n.keep*3),'b[1]'],type='l',lwd=1,col='blue') #3 chains in total


