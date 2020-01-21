#' Gillespie MCMC
#'
#' This function generates posterior samples of beta, gamma, and rho using MCMC based on Examce likelihood in subsection 4.1
#' This code refers Eq. (2.1) in Choi and Rempala (2012)
#'
#' @param indata input data set. SIR trajectory data set
#' @param nrepeat number of iteration of MCMC
#' @param prior.a hyper shape parameter of gamma prior for beta, gamma
#' @param prior.b hyper rate parameter of gamma prior for beta, gamma
#' @return returning posterior sample of beta and gamma
#' @export
Gillespie.Likelihood.MCMC <- function(indata, nrepeat, prior.a, prior.b) {
  emp.sir <- indata
  n <- emp.sir$S[1]
  time.diff <- diff(emp.sir$time)
  low.S <- emp.sir$S[-length(emp.sir$S)]
  up.S <- emp.sir$S[-1]
  low.I <- emp.sir$I[-length(emp.sir$I)]
  up.I <- emp.sir$I[-1]

  int.SI <- sum(0.5 * time.diff * (low.S * low.I + up.S * up.I))
  int.I <- sum(0.5 * time.diff * (low.I + up.I))
  zI <- (n - emp.sir$S[length(emp.sir$S)])
  zR <- emp.sir$R[length(emp.sir$R)]

  hat.rho1 <- emp.sir$I[1]/emp.sir$S[1]

  theta <- matrix(0, nrow = nrepeat, ncol = 2)
  theta[, 1] <- rgamma(nrepeat, n * zI + prior.a[1], int.SI + prior.b[1])
  theta[, 2] <- rgamma(nrepeat, zR + prior.a[2], int.I + prior.b[2])
  return(theta)
}
