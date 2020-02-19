#' SDS MCMC
#'
#' These functions are to draw posterior samples using MCMC for SDS likelihood in subsection 4.2 and Algorithm 5.1
#' This function solves the ODE in Eq. (2.7) and step 3 in Algorithm 5.1
#'
#' @param inata descrete time of an individual. This values are came from Sellke epidemic data
#' @param Tmax cutoff time of epidemic
#' @param dt time increment of ODE
#' @param beta present values of beta at the present step at MCMC
#' @param gamma present values of gamma at the present step at MCMC
#' @param ic initial values of ODE set to c(1.0, rho, 0.0), rho set the present values of beta and gamma at the present step at MCMC
#' @export
SIR.ODE <- function(indata = ti, Tmax = T.max, dt = 0.01, beta, gamma, ic) {
  p <- length(ic)
  n <- Tmax/dt
  xmat <- matrix(0, ncol = (p + 1), nrow = (n + 1))
  x <- ic
  xmat[1, ] <- c(0, x)
  jj <- 1
  for (i in 2:(n + 1)) {
    x <- x + c(-beta * x[1] * x[2], beta * x[1] * x[2] - gamma * x[2],
               gamma * x[2]) * dt
    xmat[i, 2:(p + 1)] <- x
    xmat[i, 1] <- xmat[(i - 1), 1] + dt
  }
  k <- length(indata)
  SI_ti <- matrix(0, nrow = k, ncol = 2)
  for (i in 1:k) {
    for (j in jj:nrow(xmat)) {
      if (indata[i] <= xmat[j, 1]) {
        SI_ti[i, 1:2] <- xmat[j, 2:3]
        jj <- j
        break
      }
    }
  }
  return(SI_ti)
}


#' Likelihood function for Eq. (4.3)
#'
#' Likelihood function for Eq. (4.3)
#'
#' @param SI_ti return values from SIR.ODE() function.
#' @param p.m values of beta, gamma, rho
#' @param delta duration of infectious period
#' @param n.num number of susecptible
#' @param nz.num number of removed among initially susceptible individual returning likelihood
#' @export
llikelihood <- function(SI_ti, p.m, delta.t = delta, n.num = n, nz.num = nz) {
  k <- nrow(SI_ti)
  n <- n.num
  delta <- delta.t
  lik.gamma <- nz.num * log(p.m[2]) - p.m[2] * sum(delta)
  lik <- sum(log(SI_ti[, 1])) + sum(log(SI_ti[, 2])) + k * log(p.m[1]) +
    lik.gamma + (n - k) * log(SI_ti[k, 1])
  return(lik)
}


#' SDS Likelihood MCMc
#'
#' This function generates posterior samples of beta, gamma, and rho using MCMC based on SDS likelihood in subsection 4.2
#' This function includes an additional step for generating the initial number of susceptible n. using NegativeBinomial
#' This function uses RAM method via adapt_S() function from ramcmc R-package
#'
#' @param data input data set, must be a form of Sellke epidemic (two columns. first is infection time, second is recovery time)
#' @param Tmax cutoff time of epidemic
#' @param fitn if T the function estimate initial number of susceptible, if F the function does not estimate N and calculate from data, so inputted data must cover a total population
#' @param nrepeat number of iteration of MCMC
#' @param prior.a hyper shape parameter of gamma prior for beta, gamma, rho, hyper parameter lambda
#' @param prior.b hyper rate parameter of gamma prior for beta, gamma, rho, hyper parameter lambda
#' @param ic initial value of beta, gamma, rho and N
#' @return returning posterior samples of beta, gamma, rho (with fitn=F) or beta, gamma, rho, and n (with fitn=T)
#' @export
SDS.MCMC <- function(data, Tmax, fitn = T, nrepeat = 1000,
                     prior.a, prior.b, ic = c(k1, k2, k3, k4)) {
  if(fitn == T) {
    n <- ic[4]
    theta <- matrix(0, nrow = nrepeat, ncol = 4)
  }else{
    n <- length(which(data[, 1] != 0))
    theta <- matrix(0, nrow = nrepeat, ncol = 3)
  }
  T.max <- Tmax
  delta <- c(subset((data[, 2] - data[, 1]), ((data[, 1] < T.max) & (data[, 2] < T.max))),
             subset((T.max - data[, 1])    , ((data[, 1] < T.max) & (data[, 2] >= T.max))))
  nz <- length(subset((data[, 2] - data[, 1]), ((data[, 1] < T.max) & (data[, 2] < T.max))))
  ti <- subset(data[, 1], ((data[, 1] > 0) & (data[, 1] < T.max)))
  cnt_obs = length(ti)
  parm.m <- ic[1:3]
  parm.star <- parm.m

  count <- 0;
  S <- diag(3)
  for (rep in 1:nrepeat) {
    repeat {
      u <- mvrnorm(1, c(0, 0, 0), diag(c(1, 1, 1)))
      parm.star <- parm.m + S %*% u
      if ((min(parm.star) > 0) & (parm.star[3] < 1)) {
        tau <- uniroot(function(x) 1 - x - exp(-parm.star[1]/parm.star[2] *
                                                 (x + parm.star[3])), c(0, 1))$root
        if ((tau > 0) & (tau < 1))
          break
      }
    }
    sir.m <- SIR.ODE(indata = ti, Tmax = T.max, beta = parm.m[1], gamma = parm.m[2],
                     ic = c(1, parm.m[3], 0))
    sir.star <- SIR.ODE(indata = ti, Tmax = T.max, beta = parm.star[1],
                        gamma = parm.star[2], ic = c(1, parm.star[3], 0))
    l.lik.m <- llikelihood(SI_ti = sir.m, p.m = parm.m, n.num = n,
                           delta.t = delta, nz.num = nz)
    l.lik.star <- llikelihood(SI_ti = sir.star, p.m = parm.star, n.num = n,
                              delta.t = delta, nz.num = nz)
    alpha <- exp(l.lik.star - l.lik.m
                 + dgamma(parm.star[1], prior.a[1], prior.b[1], log = T)
                 - dgamma(parm.m[1], prior.a[1], prior.b[1], log = T)
                 + dgamma(parm.star[2], prior.a[2], prior.b[2], log = T)
                 - dgamma(parm.m[2], prior.a[2], prior.b[2], log = T)
                 + dgamma(parm.star[3], prior.a[3], prior.b[3], log = T)
                 - dgamma(parm.m[3], prior.a[3], prior.b[3],   log = T))
    alpha <- min(alpha, 1)
    if (!is.nan(alpha) && runif(1) < alpha) {
      parm.m <- parm.star
      count <- count + 1
    }
    S <- ramcmc::adapt_S(S, u, alpha, rep, gamma = min(1, (3 * rep)^(-2/3)))
    theta[rep,1:3] <- parm.m

    if (fitn ==T){
      # updating tau
      x=vector('numeric',50)
      x[1]=1
      for (i in 1:(50-1)) {x[i+1]=1-exp(-parm.m[1]*(x[i]+parm.m[3])/parm.m[2])}
      tau = x[50]
      n = cnt_obs + rnbinom(1,cnt_obs,prob=tau)
      theta[rep,4] = n
    }
    if (rep%%1000 == 0) cat("0")
  }
  cat("Acceptance ratio of parameters: ",count/nrepeat,"\n")
  if(fitn ==T){
    colnames(theta) = c("beta","gamma","rho", "N")
  }else{
    colnames(theta) = c("beta","gamma","rho")
  }
  return(theta)
}


#' MCMC simulation results
#'
#' This Plot summrizes MCMC simulation results.
#'
#' @param theta inputted data set. It contains MCMC simulation results.
#' @param fitn if T the function estimate initial number of susceptible, if F the function does not estimate N and calculate from data, so inputted data must cover a total population
#' @param burn burning period
#' @param thin tinning period
#' @param nrepeat number of posterior sample
#' @return returning plots that summrizes MCMC simulation results
#' @export
result <- function(data, fitn = T, burn1=burn, nrepeat1=nrepeat, thin1=thin){
  sel = burn1+seq(0,by=thin1,length.out = (nrepeat1))
  result = data[sel,];
  stat = rbind(apply(result,2,min),
               apply(result,2,quantile,probs=0.25),
               apply(result,2,median),
               apply(result,2,mean),
               apply(result,2,quantile,probs=0.75),
               apply(result,2, max),
               apply(result,2, sd),
               apply(result,2, sd)/apply(result,2, mean))
  if (fitn==T){
    colnames(stat) = c("beta","gamma","rho", "N")
    rownames(stat) = c("Min.","1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "std", "CV")
  }else{
    colnames(stat) = c("beta","gamma","rho")
    rownames(stat) = c("Min.","1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "std", "CV")
  }
  print(stat)

  if (fitn==T){
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat = matrix(c(1, 2, 3, 4), 2, 2, byrow = T)
    layout(mat)
    plot(data[,1], type = "l", main = "beta")
    plot(data[,2], type = "l", main = "gamma")
    plot(data[,3], type = "l", main = "rho")
    plot(data[,4], type = "l", main = "N")
  }else{
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
    layout(mat)
    plot(data[,1], type = "l", main = "beta")
    plot(data[,2], type = "l", main = "gamma")
    plot(data[,3], type = "l", main = "rho")
  }
  par(mfrow=c(1,1))
}
