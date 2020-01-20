#' Sellke to trajectory
#'
#' This function converts epidemic data wih infection and removed time to SIR trajectory data with S(t), I(t), and R(t) at discrete time t
#'
#' @param sellke.data input data set, must be a form of sellke empidemic
#' @param Tmax cutoff time of epidemic
#' @return returns are SIR trajectory data set
#' @export
Sellke.to.trajectory <- function(sellke.data, Tmax) {
  mattime <- sellke.data
  n <- length(which(mattime[, 1] != 0))
  T.I <- subset(cbind(sort(mattime[, 1]), 1:length(mattime[, 1])),
                cbind(mattime[, 1] >= 0))
  colnames(T.I) <- c("time", "I")
  T.I <- as.data.frame(T.I)
  T.R <- subset(cbind(sort(mattime[, 2]), 1:length(mattime[, 2])),
                cbind(mattime[, 2] >= 0))
  colnames(T.R) <- c("time", "R")
  T.R <- as.data.frame(T.R)

  all <- merge(x = T.I, y = T.R, by = "time", all = T)
  all.1 <- subset(all, all[, 1] == 0)
  all.2 <- subset(all, all[, 1] != 0)
  all.1 <- all.1[nrow(all.1), ]
  all <- rbind(all.1, all.2)
  rm(all.1, all.2)
  S <- rep(n, nrow(all))
  all <- cbind(all, S)
  for (i in 2:nrow(all)) {
    if (is.na(all$I[i]))
      all$I[i] <- all$I[i - 1]
    if (is.na(all$R[i]))
      all$R[i] <- all$R[i - 1]
    if (is.na(all$R[i]))
      all$R[i] <- 0
    all$S[i] <- all$S[1] - all$I[i] + all$I[1]
  }
  all$R[1] <- 0
  all$I <- all$I - all$R
  emp.sir <- subset(all, all$time < Tmax)
  emp.sir <- cbind(emp.sir, 1:nrow(emp.sir))
  time.diff <- diff(emp.sir$time)
  for (i in 1:(nrow(emp.sir) - 1)) {
    if ((time.diff[i] == 0) & emp.sir[i, 1] == 0)
      emp.sir[i, 5] <- 0
  }
  emp.sir <- subset(emp.sir, emp.sir[, 5] != 0)
  emp.sir <- emp.sir[, 1:4]
  return(emp.sir)
}
