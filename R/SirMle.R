#' SIR MLE
#'
#' This function calculates MLE using SIR empidemic data
#'
#' @param indata input data set. SIR trajectory data set
#' @return returning MLE of beta, gamma, rho of subsection 4.1
#' @export
SIR.MLE <- function(indata) {
  emp.sir <- indata
  n <- emp.sir$S[1]
  time.diff <- diff(emp.sir$time)
  low.S <- emp.sir$S[-length(emp.sir$S)]
  up.S <- emp.sir$S[-1]
  low.I <- emp.sir$I[-length(emp.sir$I)]
  up.I <- emp.sir$I[-1]

  int.SI <- sum(0.5 * time.diff * (low.S * low.I + up.S * up.I))
  int.I <- sum(0.5 * time.diff * (low.I + up.I))

  hat.beta1 <- n * (n - emp.sir$S[length(emp.sir$S)])/int.SI
  hat.gamma1 <- emp.sir$R[length(emp.sir$R)]/int.I
  hat.rho1 <- emp.sir[1, 2]/n
  return(c(hat.beta1, hat.gamma1, hat.rho1))
}
