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

# Third example data: simulation data. The data is simulated by Sellke construction.
# This example has two MCMC simulations using synthetic data by Sellke construction. 
# Running the other methods on the paper except SDS approach

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

#converting Sellke epidemic data to SIR trajectory
emp.sir <- Sellke.to.trajectory(pop.data, Tmax = T.max)

# MLE using Mehtod 1
MLE <- SIR.MLE(emp.sir)
print(MLE)
# MCMC using Method 2
gillespie <- Gillespie.Likelihood.MCMC(emp.sir, nrepeat=sim.num, prior.a=c(0.001,0.001), prior.b=c(0.001,0.001)) 

summary(gillespie)
plot(gillespie[,1],type="l")
plot(gillespie[,2],type="l")


#MCMC using Method 4
gaussian <- Gaussian.MCMC(emp.sir, nrepeat=sim.num, ic = c(1, 1, 0.1), tun=c(0.1,0.1,1), prior.a=c(0.001,0.001,0.05*0.1)
                          , prior.b=c(0.001,0.001,0.1), Tmax = T.max)
result(data = gaussian, fitn = F)

