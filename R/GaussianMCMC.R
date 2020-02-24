#' Gaussian MCMC
#'
#' This function generates posterior samples of beta, gamma, and rho using MCMC based on Gaussian likelihood in subsection A.2
#'
#' @param data input data set. SIR trajectory data set
#' @param nrepeat number of iteration of MCMC
#' @param ic initial value of beta, gamma, and rho
#' @param tun tunning constant for proposal distribution of beta, gamma, rho
#' @param prior.a hyper shape parameter of gamma prior for beta, gamma, rho
#' @param prior.b hyper rate parameter of gamma prior for beta, gamma, rho
#' @param Tmax cutoff time of epidemic
#' @return returning posterior samples of beta, gamma, rho
#' @export
Gaussian.MCMC <- function(data, nrepeat, ic = c(k1, k2, k3), tun, prior.a, prior.b,
                          Tmax = 100) {
  T.max <- Tmax
  emp.sir <- data
  n <- emp.sir$S[1]
  time.diff <- diff(emp.sir$time)
  low.S <- emp.sir$S[-length(emp.sir$S)]
  up.S <- emp.sir$S[-1]
  low.I <- emp.sir$I[-length(emp.sir$I)]
  up.I <- emp.sir$I[-1]

  int.SI <- sum(0.5 * time.diff * (low.S * low.I + up.S * up.I))
  int.I <- sum(0.5 * time.diff * (low.I + up.I))
  Z.I <- (n - emp.sir$S[length(emp.sir$S)])
  Z.R <- emp.sir$R[length(emp.sir$R)]
  mu1.beta <- Z.I/(int.SI/n)
  mu0.beta <- prior.a[1]/prior.b[1]
  c0.beta <- 1/(prior.a[1]/(prior.b[1]**2))
  mu1.gamma <- Z.R/(int.I)
  mu0.gamma <- prior.a[2]/prior.b[2]
  c0.gamma <- 1/(prior.a[2]/(prior.b[2]**2))
  parm.m <- ic

  theta <- matrix(0, nrow = nrepeat, ncol = 3)
  count <- 0
  for (rep in 1:nrepeat) {
    repeat {
      a.star <- rnorm(1, mean = parm.m[3], tun[3])
      if ((a.star > 0) & (a.star < 1))
        break
    }

    sir.m <- SIR.ODE.full(Tmax = T.max, beta = parm.m[1], gamma = parm.m[2],
                          ic = c(1, parm.m[3], 0), dt = 0.1)
    sir.star <- SIR.ODE.full(Tmax = T.max, beta = parm.m[1], gamma = parm.m[2],
                             ic = c(1, a.star, 0), dt = 0.1)
    row <- nrow(sir.m)
    l.lik.m <- (dnorm(Z.I, parm.m[1] * int.SI/n, sqrt(n * (1 - sir.m[row, 2])), log = T)
                    + dnorm(Z.R, parm.m[2] * int.I, sqrt(n * (1 +
                      parm.m[3] - sir.m[row, 2] - sir.m[row, 3])), log = T))
    l.lik.star <- (dnorm(Z.I, parm.m[1] * int.SI/n, sqrt(n * (1 - sir.star[row,2]))
                   , log = T) + dnorm(Z.R, parm.m[2] * int.I, sqrt(n * (1 +
                   a.star - sir.star[row, 2] - sir.star[row, 3])), log = T))
    alpha <- exp(l.lik.star + dgamma(a.star, prior.a[3], prior.b[3], log = T)
                 + pnorm(a.star, 0, tun[3], log.p = T)
                 - l.lik.m - dgamma(parm.m[3], prior.a[3], prior.b[3], log = T)
                 - pnorm(parm.m[3], 0, tun[3], log.p = T))
    if (!is.nan(alpha) && runif(1) < alpha) {
      parm.m[3] <- a.star
      sir.m <- SIR.ODE.full(Tmax = T.max, beta = parm.m[1], gamma = parm.m[2],
                            ic = c(1, parm.m[3], 0), dt = 0.1)
      count <- count + 1
    }
    c1.beta <- 1/((n * (1 - sir.m[row, 2]))/((int.SI/n)^2))
    c1.gamma <- 1/((n * (1 + parm.m[3] - sir.m[row, 2] - sir.m[row,
                   3]))/((int.I)^2))
    sig.beta <- 1/(c0.beta + c1.beta)
    sig.gamma <- 1/(c0.gamma + c1.gamma)
    mu.beta <- sig.beta * (c1.beta * mu1.beta + c0.beta * mu0.beta)
    mu.gamma <- sig.gamma * (c1.gamma * mu1.gamma + c0.gamma * mu0.gamma)

    repeat {
      parm.m[1] <- rnorm(1, mean = mu.beta, sd = sqrt(sig.beta))
      if (parm.m[1] > 0)
        break
    }
    repeat {
      parm.m[2] <- rnorm(1, mean = mu.gamma, sd = sqrt(sig.gamma))
      if (parm.m[2] > 0)
        break
    }
    theta[rep, ] <- parm.m
    if (rep%%100 == 0)
      cat("0")  #cat(parm.m,'\n')
  }
  print(count/nrepeat)
  return(theta)
}


#' This function solves the ODE in Eq. (2.7)
#'
#' This function solves the ODE in Eq. (2.7)
#'
#' @param Tmax cutoff time of epidemic
#' @param dt time increment of ODE
#' @param beta present values of beta at the present step at MCMC
#' @param gamma present values of gamma at the present step at MCMC
#' @param ic initial values of ODE set to c(1.0, rho, 0.0), rho set the present values of beta and gamma at the present step at MCMC
#' @return
#' @export
SIR.ODE.full <- function(Tmax, dt = 0.01, beta, gamma, ic) {
  p = length(ic)
  n = Tmax/dt
  xmat = matrix(0, ncol = (p + 1), nrow = (n + 1))
  x = ic
  xmat[1, ] = c(0, x)
  for (i in 2:(n + 1)) {
    x = x + c(-beta * x[1] * x[2], beta * x[1] * x[2] - gamma * x[2],
              gamma * x[2]) * dt
    xmat[i, 2:(p + 1)] = x
    xmat[i, 1] = xmat[(i - 1), 1] + dt
  }
  return(xmat)
}
