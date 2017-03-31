# rm(list=ls())
# require(R2jags)
library(jagsUI)
library(MCMCpack) # rwish function
library(plyr)
library(data.table)
library(car)
library(Metrics)

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
cpe_dat <- cpe_dat[, c("dowlknum","Year","walleye","yellow.perch","northern.pike","DOW"), with=FALSE]
# log-transfrom CPE
cpe_dat[, log_yp_cpe:= log(yellow.perch+0.1)]
cpe_dat[, log_wae_cpe:= log(walleye+0.1)]
cpe_dat[, log_np_cpe:= log(northern.pike+0.1)]

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
cpe_gdd[, z_np_cpe:= as.numeric(scale(log_np_cpe))]
cpe_gdd[, z_gdd:= as.numeric(scale(gdd))]
cpe_gdd[, z_gdd1:= as.numeric(scale(gdd.lag1))]
cpe_gdd[, z_gdd2:= as.numeric(scale(gdd.lag2))]
cpe_gdd[, z_gdd3:= as.numeric(scale(gdd.lag3))]
cpe_gdd[, z_gdd4:= as.numeric(scale(gdd.lag4))]
cpe_gdd[, z_gdd5:= as.numeric(scale(gdd.lag5))]
cpe_gdd[, z_gddma:= as.numeric(scale(gdd.ma.5))]

head(cpe_gdd)
dim(cpe_gdd)
# Correlation among gdd predictors and wae
gdd_cor <- cpe_gdd[,c(17:25)]
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

# Create lake-level predictors of WAE and NP CPE
wae_lake <- as.numeric(by(cpe_gdd$log_wae_cpe, cpe_gdd$DOW, mean)) 
length(wae_lake)
z_wae_lake <- as.numeric(scale(wae_lake))

np_lake <- as.numeric(by(cpe_gdd$log_np_cpe, cpe_gdd$DOW, mean)) 
length(np_lake)
z_np_lake <- as.numeric(scale(np_lake))

# hist(lake_dat_sub$alkalinity)
#################################################################
########## BUGS CODE ############################################
#################################################################
# Define the model in the BUGS language and write a text file
sink("model.txt")
cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]] + beta1[group[i]] * x1[i] + beta2[group[i]] * x2[i] + beta3[group[i]] * x3[i]        
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta1[j] <- BB[j,2]
    beta2[j] <- BB[j,3]
    beta3[j] <- BB[j,4]

    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivriate normal
    
    BB.hat[j,1] <- mu.alpha + s1[1] * z1[j] + s1[2] * z2[j] + s1[3] * z3[j] + s1[4] * z4[j] 
      + s1[5] * z5[j] + s1[6] * z6[j] 
      BB.hat[j,2] <- mu.beta1 + s2[1] * z1[j] + s2[2] * z2[j] + s2[3] * z3[j] + s2[4] * z4[j] 
      + s2[5] * z5[j] + s2[6] * z6[j] 
      BB.hat[j,3] <- mu.beta2 + s3[1] * z1[j] + s3[2] * z2[j] + s3[3] * z3[j] + s3[4] * z4[j] 
      + s3[5] * z5[j] + s3[6] * z6[j] 
      BB.hat[j,4] <- mu.beta3 + s4[1] * z1[j] + s4[2] * z2[j] + s4[3] * z3[j] + s4[4] * z4[j] 
      + s4[5] * z5[j] + s4[6] * z6[j] 
    }
    
    
    # Bayesian LASSO -  a Laplace (double exponential) prior
    # Level-2 predictors
    for(k in 1:6){
    s1[k] ~ ddexp(0, lambda1)
    }
    
    for(k in 1:6){
    s2[k] ~ ddexp(0, lambda2)
    }

    for(k in 1:6){
    s3[k] ~ ddexp(0, lambda3)
    }

    for(k in 1:6){
    s4[k] ~ ddexp(0, lambda4)
    }
    
    # Hyper-prior on lambda
    lambda1 ~ dgamma(0.1,0.1)
    lambda2 ~ dgamma(0.1,0.1)
    lambda3 ~ dgamma(0.1,0.1)
    lambda4 ~ dgamma(0.1,0.1)
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    tau <- pow(sigma,-2) # precision
    sigma2 <- pow(sigma,2)
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta1 ~ dnorm(0, 0.0001)
    mu.beta2 ~ dnorm(0, 0.0001)
    mu.beta3 ~ dnorm(0, 0.0001)
    
    
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()



head(cpe_gdd)
# Number of sites
J <-length(unique(cpe_gdd$DOW))

# Site indicator
cpe_gdd[, G := as.numeric(as.factor(as.numeric(DOW)))]

# Number of varying parameters
K <- 4

# Create identity matrix for Wishart dist'n
#!!!!!!!Number of parameters to estimate (K)

W <- diag(K)

# Load data
data <- list(y = cpe_gdd$log_yp_cpe, group = cpe_gdd$G, n = cpe_gdd[,.N], J = J,
             x1 = cpe_gdd$z_gdd1, x2 = cpe_gdd$z_wae_cpe, x3 = cpe_gdd$log_np_cpe,
             z1 = lake_dat_sub$z_area, z2 = lake_dat_sub$z_littoral, z3 = lake_dat_sub$z_depth, z4 = lake_dat_sub$z_gddMean,
             z5 = z_wae_lake, z6 = z_np_lake,
             K=K,W=W)




# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1), sigma=runif(1), sigma.alpha=runif(1) )
}


# Parameters monitored
parameters <- c("mu.alpha","sigma", "sigma.alpha","b","lambda1","lambda2","lambda3","lambda4",
                "s1","s2","s3","s4","mu.beta1","mu.beta2","mu.beta3")


# MCMC settings
ni <- 110000
nt <- 2
nb <- 80000
nc <- 3


start.time = Sys.time()         # Set timer 
# Call JAGS from R 

out <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time

# Find which parameters, if any, have Rhat > 1.1
which(out$summary[, c("Rhat")] > 1.1)


# Summarize posteriors
print(out, dig = 3)
# traceplot(out)
# outExp2 <- out$summary
# write.csv(outExp2, "CSI_ModelSummary.csv", row.names = T)
mcmcOut <- out$sims.list
saveRDS(mcmcOut, file="CSI_mcmc_out.rds")

# str(out)


sEsts1 <- matrix(NA, nrow=3,ncol=6)
for(i in 1:6){
  sEsts1[1,i] <- mean(out$sims.list$s1[,i])
  sEsts1[2:3,i] <- quantile(out$sims.list$s1[,i],c(0.025,0.975))
}
sEsts1

sEsts2 <- matrix(NA, nrow=3,ncol=6)
for(i in 1:6){
  sEsts2[1,i] <- mean(out$sims.list$s2[,i])
  sEsts2[2:3,i] <- quantile(out$sims.list$s2[,i],c(0.025,0.975))
}
sEsts2

sEsts3 <- matrix(NA, nrow=3,ncol=6)
for(i in 1:6){
  sEsts3[1,i] <- mean(out$sims.list$s3[,i])
  sEsts3[2:3,i] <- quantile(out$sims.list$s3[,i],c(0.025,0.975))
}
sEsts3

sEsts4 <- matrix(NA, nrow=3,ncol=6)
for(i in 1:6){
  sEsts4[1,i] <- mean(out$sims.list$s4[,i])
  sEsts4[2:3,i] <- quantile(out$sims.list$s4[,i],c(0.025,0.975))
}
sEsts4


Ests <- cbind(sEsts1, sEsts2, sEsts3, sEsts4)

#####################################################
########### PLOT ####################################
#####################################################


covariates <- c("Lake area", "Littoral area", "Depth", "Mean GDD", "Mean WAE", "Mean NP")

res <- 6
name_figure <- "CSI_MN_YP_effects.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

#nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
# nf <- layout(matrix(c(1:2), nrow = 1, ncol=2,,byrow=TRUE), widths = c(0.7, 0.3))
nf <- layout(matrix(c(1:4), nrow = 2, ncol=2,byrow=TRUE))
layout.show(nf)
#par(mar=c(1,1,1,1), oma=c(2,2,1,1) )
par(mar = c(1, 3.7, 0, 0) + 0.1,oma=c(2,1.5,0,0))
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

size.labels = 1
size.text = 1

x.label <- 'Estimated effect'
y.label <- 'Covariate'

# Posterior means and CIs for all parameters
Plot.data <- Ests[,1:6]

rows <- 1:dim(Plot.data)[2]

Plot.color <- as.numeric(Plot.data[2,] * Plot.data[3,] > 0 )
colorP <- rep("black", length(rows))
colorP[Plot.color > 0] <- "blue"

# pchP <- rep(16, length(covariates))
# pchP[10:15] <- 15

pchP <- 15

plotting.region <- range(Ests)



### axis label options
spc <- 0.26
lab <- 1:63
cex <- 0.5
adj <- 0
###
plot(c(plotting.region[1], plotting.region[2]), c(0.5,length(rows)+0.5), 
     axes=F, xlab='',ylab='',type='n')
axis(side=1,cex.axis=size.text, mgp=c(0,0.1,0),tck= -0.01, labels=F) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:length(rows)),labels=F,tck= -0.01)

text(par("usr")[1] - spc,1:length(rows),srt = 0, adj =adj,
     labels = covariates, xpd = TRUE,cex=0.7)

# 95% CIs for censored analysis
segments(x0=Plot.data[2,], x1=Plot.data[3,],
         y0=1:length(rows), y1=1:length(rows), col=colorP,lwd=1)

# segments(x0=Plot.data[,4], x1=Plot.data[,6],
#          y0=1:length(rows), y1=1:length(rows), lwd=4, col="blue")
## Estiamtes from censored model
points(Plot.data[1,], 1:length(rows), col=colorP ,cex=1, pch=pchP)

abline(v=0, col="gray")
# abline(v=out2$mean$mu.ave, col="gray")

# Add x- and y-axis lables
# mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T, adj=0.6)
# mtext(y.label, line = 0.4, side = 2, cex = size.text, outer=T)
text(-0.4,6.4,"Lake mean YP CEP", cex=0.8)


# abline(h=0)
box()


############## PLOT 2
# Posterior means and CIs for all parameters
Plot.data <- Ests[,7:12]

par(mar = c(1, 1, 0, 0) )

rows <- 1:dim(Plot.data)[2]

Plot.color <- as.numeric(Plot.data[2,] * Plot.data[3,] > 0 )
colorP <- rep("black", length(rows))
colorP[Plot.color > 0] <- "blue"

# pchP <- rep(16, length(covariates))
# pchP[10:15] <- 15



### axis label options
spc <- 0.26
lab <- 1:63
cex <- 0.5
adj <- 0
###
plot(c(plotting.region[1], plotting.region[2]), c(0.5,length(rows)+0.5), 
     axes=F, xlab='',ylab='',type='n')
axis(side=1,cex.axis=size.text, mgp=c(0,0.1,0),tck= -0.01, labels=F) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:length(rows)),labels=F,tck= -0.01)

# text(par("usr")[1] - spc,1:length(rows),srt = 0, adj =adj,
#      labels = covariates, xpd = TRUE,cex=0.7)

# 95% CIs for censored analysis
segments(x0=Plot.data[2,], x1=Plot.data[3,],
         y0=1:length(rows), y1=1:length(rows), col=colorP,lwd=1)

# segments(x0=Plot.data[,4], x1=Plot.data[,6],
#          y0=1:length(rows), y1=1:length(rows), lwd=4, col="blue")
## Estiamtes from censored model
points(Plot.data[1,], 1:length(rows), col=colorP ,cex=1, pch=pchP)

abline(v=0, col="gray")
# abline(v=out2$mean$mu.ave, col="gray")

# Add x- and y-axis lables
# mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T, adj=0.6)
# mtext(y.label, line = 0.4, side = 2, cex = size.text, outer=T)

text(-0.32,6.4,"Effect of GDD (lag 1) on YP CEP", cex=0.8)


# abline(h=0)
box()


############## PLOT 3
# Posterior means and CIs for all parameters
Plot.data <- Ests[,13:18]

par(mar = c(1, 3.7, 0, 0) )

rows <- 1:dim(Plot.data)[2]

Plot.color <- as.numeric(Plot.data[2,] * Plot.data[3,] > 0 )
colorP <- rep("black", length(rows))
colorP[Plot.color > 0] <- "blue"

# pchP <- rep(16, length(covariates))
# pchP[10:15] <- 15




### axis label options
spc <- 0.255
lab <- 1:63
cex <- 0.5
adj <- 0
###
plot(c(plotting.region[1], plotting.region[2]), c(0.5,length(rows)+0.5), 
     axes=F, xlab='',ylab='',type='n')
axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:length(rows)),labels=F,tck= -0.01)

text(par("usr")[1] - spc,1:length(rows),srt = 0, adj =adj,
     labels = covariates, xpd = TRUE,cex=0.7)

# 95% CIs for censored analysis
segments(x0=Plot.data[2,], x1=Plot.data[3,],
         y0=1:length(rows), y1=1:length(rows), col=colorP,lwd=1)

# segments(x0=Plot.data[,4], x1=Plot.data[,6],
#          y0=1:length(rows), y1=1:length(rows), lwd=4, col="blue")
## Estiamtes from censored model
points(Plot.data[1,], 1:length(rows), col=colorP ,cex=1, pch=pchP)

abline(v=0, col="gray")
# abline(v=out2$mean$mu.ave, col="gray")

# Add x- and y-axis lables
# mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T, adj=0.6)
# mtext(y.label, line = 0.4, side = 2, cex = size.text, outer=T)

text(-0.34,6.4,"Effect of WAE on YP CEP", cex=0.8)


# abline(h=0)
box()


############## PLOT 4
# Posterior means and CIs for all parameters
Plot.data <- Ests[,19:24]

par(mar = c(1, 1, 0, 0) + 0.1)


rows <- 1:dim(Plot.data)[2]

Plot.color <- as.numeric(Plot.data[2,] * Plot.data[3,] > 0 )
colorP <- rep("black", length(rows))
colorP[Plot.color > 0] <- "blue"

# pchP <- rep(16, length(covariates))
# pchP[10:15] <- 15


# plotting.region <- range(Plot.data)

### axis label options
spc <- 0.001
lab <- 1:63
cex <- 0.5
adj <- 0
###
plot(c(plotting.region[1], plotting.region[2]), c(0.5,length(rows)+0.5), 
     axes=F, xlab='',ylab='',type='n')
axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:length(rows)),labels=F,tck= -0.01)

# text(par("usr")[1] - spc,1:length(rows),srt = 0, adj =adj,
#      labels = covariates, xpd = TRUE,cex=0.7)

# 95% CIs for censored analysis
segments(x0=Plot.data[2,], x1=Plot.data[3,],
         y0=1:length(rows), y1=1:length(rows), col=colorP,lwd=1)

# segments(x0=Plot.data[,4], x1=Plot.data[,6],
#          y0=1:length(rows), y1=1:length(rows), lwd=4, col="blue")
## Estiamtes from censored model
points(Plot.data[1,], 1:length(rows), col=colorP ,cex=1, pch=pchP)

abline(v=0, col="gray")
# abline(v=out2$mean$mu.ave, col="gray")

# Add x- and y-axis lables
mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T, adj=0.5)
mtext(y.label, line = 0.4, side = 2, cex = size.text, outer=T)

text(-0.4,6.4,"Effect of NP on YP CEP", cex=0.8)

# abline(h=0)
box()

par(def.par)
dev.off()



