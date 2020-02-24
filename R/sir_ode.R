#'ode.final.size
#'
#'euler.prob
#'rho=.038
#'beta=0.2160
#'gamma=0.1928
#'
#' @param rho
#' @param gamma
#' @param lambda
#' @param k
#' @return
#' @export
ode.final.size=function(rho=0.01, gamma=1,lambda=2.9,k=20)
{ x=vector('numeric',k)
x[1]=1
for (i in 1:(k-1)) {x[i+1]=1-exp(-lambda*(x[i]+rho)/gamma)}
return(x) }
#ode.final.size(0.02,.25,.23)


#'ode.final.size.cnt
#'
#'euler.prob
#'rho=.038
#'beta=0.2160
#'gamma=0.1928
#'
#' @param x0
#' @param y0
#' @param gamma
#' @param lambda
#' @param k
#' @return
#' @export
ode.final.size.cnt=function(x0=100, y0=1, gamma=1,lambda=1.9,k=20)
{x=vector('numeric',k)
n=x0+y0
x[1]=1 #sample(1:round(x0))[1]
for (i in 1:(k-1)) {x[i+1]=1-x0*exp(-lambda*x[i]/gamma)/n}
return(x*n) }
#ode.final.size.cnt()


#'euler
#'
#'euler.prob
#'rho=.038
#'beta=0.2160
#'gamma=0.1928
#'
#' @param t
#' @param dt
#' @param fun
#' @param ic
#' @return
#' @export
euler <- function(t = 120, dt = 0.001, fun = f, ic=c(1,rho,0))
{
p <- length(ic)
 n <- t/dt
 xmat <- matrix(0,ncol=p,nrow=n)
 x <- ic
 xmat[1,] <- x
 for (i in 2:n) {
 	x <- x + fun(x)*dt
	xmat[i,] <- x
 }
ts(xmat,start=0,deltat=dt)
}


#'ode.sir
#'
#'euler.prob
#'rho=.038
#'beta=0.2160
#'gamma=0.1928
#'
#' @param x
#' @param k1
#' @param k2
#' @return
#' @export
ode.sir <- function(x, k1=beta, k2=gamma)
 {
  c( -k1*x[1]*x[2] ,
  k1*x[1]*x[2] - k2*x[2], k2*x[2])
  }
#' @example plot.ts(euler(fun=ode.sir),plot.type="si",ylab='SIR_ODE')
