#' Sellke
#'
#' This function generate synthetic epidemic data using Sellke construction in Algorithm 3.1
#'
#' @param n number of initial susceptible
#' @param beta parmeters, tranmission rate, recovery rate, initial infective ratio
#' @param gamma parmeters, tranmission rate, recovery rate, initial infective ratio
#' @param rho parmeters, tranmission rate, recovery rate, initial infective ratio
#' @param Tmax cutoff time of epidemic
#' @param dt time increment
#' @return returns are individual infection time and removed time.
#' @return number of rows are n + n*rho
#' @export
Sellke <- function(n, rho = 0.01, beta = 0.3, gamma = 0.15, Tmax = 80,
                   dt = 0.001) {
  # n=pop; rho=k3; beta=k1; gamma=k2; Tmax = T.max;dt=0.001;
  m <- round(n * rho)
  Q <- sort(rexp(n, 1))
  delta <- rexp((n + m), gamma)


  # matrix for saving infection time and removal for each individual
  mattime <- matrix(NA, nrow = (n + m), ncol = 2)
  for (i in 1:m) {
    mattime[i, 1] <- 0
    mattime[i, 2] <- 0 + delta[i]
  }

  I <- m  # initial number of infection
  R <- 0  # initial number of recover
  t <- 0
  At <- 0
  min.pre <- 0
  min.T.R <- min(subset(mattime[, 2], mattime[, 2] > min.pre))
  for (i in 1:n) {
    for (time in seq(t, Tmax, by = dt)) {
      if (time >= min.T.R) {
        I <- I - 1
        R <- R + 1
        t <- time
        min.pre <- min.T.R
        if (I == 0)
          break
        min.T.R <- min(subset(mattime[, 2], mattime[, 2] > min.pre))
      }
      if (At > Q[i]) {
        I <- I + 1
        t <- time
        mattime[m + i, 1] <- t
        mattime[m + i, 2] <- t + delta[m + i]
        min.T.R <- min(subset(mattime[, 2], mattime[, 2] > min.pre))
        break
      }
      At <- At + (beta/n) * I * dt
    }
    if (I == 0)
      break
  }
  mattime1 <- subset(mattime, is.na(mattime[, 1]) == F)
  mattime2 <- subset(mattime, is.na(mattime[, 1]) == T)
  mattime2[, 1] <- Tmax
  mattime2[, 2] <- mattime2[, 1] + rexp(nrow(mattime2), gamma)
  in.data <- rbind(mattime1, mattime2)
  return(in.data)
}

