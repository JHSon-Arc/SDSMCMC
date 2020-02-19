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
# The first example data: H1N1 pandemic data from WSU campus 
# In this example, we estimate three parameters and the initial number of susceptible 
# Input data should have two columns. The first is infection time and the second is recovery time. 


data = read.csv(file = "WSU.csv")
T.max=104; # cut-off time point 
n1=0; # initial number of susceptible 
burn = 1000; #  burning period
thin =1; # tinning period 
nrepeat = 2000; # number of posterior sample;  
sim.num= burn + nrepeat * thin; # total number of simulation 

sds  <- SDS.MCMC(data = data, Tmax=T.max, fitn = T, nrepeat = sim.num, 
                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(0.25, 0.2, 0.01, 1000))
result(data = sds, fitn = T)
