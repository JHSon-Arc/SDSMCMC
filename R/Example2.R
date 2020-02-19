# data from UW epidemic as described in H1N1 paper with ES 

library(plyr)
library(MASS)
library(smfsb)
library(ramcmc) 
#library(SDSMCMC)

source("SdsMCMC.R")
source("GaussianMCMC.R")
source("GillespieMCMC.R")
source("Sellke.R")
source("SellkeToTrajecotory.R")
source("SirMle.R")
source("sir_ode.R")

############################################################################
# Second example data: simulation data. The data is simulated by Sellke construction.
# This example has two MCMC simulations using synthetic data by Sellke construction. 
#The first is estimation including estimation n, number of susceptible. The second runs without n.

burn = 1000; #  burning period
thin =1; # tinning period 
nrepeat = 2000; # number of posterior sample;  
sim.num= burn + nrepeat * thin; # total number of simulation 

#initial parameter setting
k1 = 1.1; k2 = 0.8 ; k3 = 0.05; n=1000; T.max = 30; 
beta=k1; gamma=k2; rho=k3;
plot.ts(euler(fun=ode.sir),plot.type="si",ylab='SIR_ODE')

#generating synthetic epidemic data using Sellke construction
pop.data = Sellke(n=n, rho=k3, beta=k1, gamma=k2, Tmax = T.max)

#MCMC using Method 3. In this example, we estimate three parameters given n = 1000 
sds  <- SDS.MCMC(fitn = F, data = pop.data, Tmax=T.max, nrepeat = sim.num, 
                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(k1, k2, k3, n))
result(data = sds, fitn = F)


# MCMC using Method 3. In this example, we estimate three parameters and the initial number of susecptible 
sds  <- SDS.MCMC(fitn = T, data = pop.data, Tmax=T.max, nrepeat = sim.num, 
                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(k1, k2, k3, n))
result(data = sds, fitn = T)
