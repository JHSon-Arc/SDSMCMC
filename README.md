# SDSMCMC

Unidirectional Mass Transfer Models (MTMs) are compartmental dynamical systems with a notion of partial ordering among the compartments. One can interpret a unidirectional MTM as a Survival Dynamical System (SDS) that is fully described in terms of survival functions (as opposed to population counts in MTMs). In this repository, we provide the necessary toolkit for a statistical inference method based on the SDS interpretation of the compartmental susceptible-infected-recovered (SIR) epidemic model due to Kermack and McKendrick for the papaer on https://arxiv.org/abs/1207.3137.

We provide 7 R scripts.

Sellke.R: This function generate synthetic epidemic data using Sellke construciotn in Algorithm 3.1.

SellkeToTrajectory.R: This function converts epidemic data wih infection and removed time to SIR trajectory data with S(t), I(t), and R(t) at discrete time t

SirMle.R: This function calculates MLE using SIR empidemic data.

GillespieMCMC.R: This function generates posterior samples of beta, gamma, and rho using MCMC based on Examce likelihood in subsection 4.1.

GaussianMCMC.R: This function generates posterior samples of beta, gamma, and rho using MCMC based on Gaussian likelihood in subsection A.2.

SdsMCMC.R: These functions are to draw posterior samples using MCMC for SDS likelihood in subsection 4.2 and Algorithm 5.1.

<pre>
# Example1
<code>
# data from UW epidemic as described in H1N1 paper with ES 

library(plyr)
library(MASS)
library(smfsb)
library(ramcmc) 
library(SDSMCMC)

# source("SdsMCMC.R")
# source("GaussianMCMC.R")
# source("GillespieMCMC.R")
# source("Sellke.R")
# source("SellkeToTrajecotory.R")
# source("SirMle.R")
# source("sir_ode.R")

############################################################################
# The first example data: H1N1 pandemic data from WSU campus 
# In this example, we estimate three parameters and the initial number of susceptible 
# Input data should have two columns. The first is infection time and the second is recovery time. 


data = WSU
T.max=104; # cut-off time point 
n1=0; # initial number of susceptible 
burn = 1000; #  burning period
thin =1; # tinning period 
nrepeat = 5000; # number of posterior sample;  
sim.num= burn + nrepeat * thin; # total number of simulation 

sds  <- SDS.MCMC(data = data, Tmax=T.max, fitn = T, nrepeat = sim.num, 
                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(0.25, 0.2, 0.01, 1000))
result(data = sds, fitn = T)
</code>
</pre>

<pre>
# Example2
<code>
# data from UW epidemic as described in H1N1 paper with ES 

library(plyr)
library(MASS)
library(smfsb)
library(ramcmc) 
library(SDSMCMC)

# source("SdsMCMC.R")
# source("GaussianMCMC.R")
# source("GillespieMCMC.R")
# source("Sellke.R")
# source("SellkeToTrajecotory.R")
# source("SirMle.R")
# source("sir_ode.R")

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

#generating synthetic epidemic data using Sellke construction
pop.data = Sellke(n=n, rho=k3, beta=k1, gamma=k2, Tmax = T.max)

#MCMC using Method 3. In this example, we estimate three parameters given n = 1000 
sds  <- SDS.MCMC(data = pop.data, Tmax=T.max, fitn = F, nrepeat = sim.num, 
                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(k1, k2, k3, n))
result(data = sds, fitn = F)


# MCMC using Method 3. In this example, we estimate three parameters and the initial number of susecptible 
sds  <- SDS.MCMC(data = pop.data, Tmax=T.max, fitn = T, nrepeat = sim.num, 
                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(k1, k2, k3, n))
result(data = sds, fitn = T)
</code>
</pre>

<pre>
# Example3
<code>
# data from UW epidemic as described in H1N1 paper with ES 

library(plyr)
library(MASS)
library(smfsb)
library(ramcmc) 
library(SDSMCMC)

#source("SdsMCMC.R")
#source("GaussianMCMC.R")
#source("GillespieMCMC.R")
#source("Sellke.R")
#source("SellkeToTrajecotory.R")
#source("SirMle.R")
#source("sir_ode.R")

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
</code>
</pre>

