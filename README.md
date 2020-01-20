# SDSMCMC

Unidirectional Mass Transfer Models (MTMs) are compartmental dynamical systems with a notion of partial ordering among the compartments. One can interpret a unidirectional MTM as a Survival Dynamical System (SDS) that is fully described in terms of survival functions (as opposed to population counts in MTMs). In this repository, we provide the necessary toolkit for a statistical inference method based on the SDS interpretation of the compartmental susceptible-infected-recovered (SIR) epidemic model due to Kermack and McKendrick for the papaer on https://arxiv.org/abs/1207.3137.

We provide 7 R scripts.

Sellke.R: This function generate synthetic epidemic data using Sellke construciotn in Algorithm 3.1.

SellkeToTrajectory.R: This function converts epidemic data wih infection and removed time to SIR trajectory data with S(t), I(t), and R(t) at discrete time t

SirMle.R: This function calculates MLE using SIR empidemic data.

GillespieMCMC.R: This function generates posterior samples of beta, gamma, and rho using MCMC based on Examce likelihood in subsection 4.1.

GaussianMCMC.R: This function generates posterior samples of beta, gamma, and rho using MCMC based on Gaussian likelihood in subsection A.2.

SdsMCMC.R: These functions are to draw posterior samples using MCMC for SDS likelihood in subsection 4.2 and Algorithm 5.1.

Example.R: This file provides examples using R codes to generate synthetic epidemic data and to estimate parameters accoriding to estimation methods described in the paper.

WS_data_MCMC_NB.R: R code of MCMC simulation for H1N1 pandemic data at Washington State University data1.csv: Daily count data of H1N1 pandemic at Washington State University
