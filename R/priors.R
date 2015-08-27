#' Creates expert prior from ellicitation
#' @param m Median for 1\% value at risk
#' @param q 90\% quantile for median m
#' @param mDiff Median for difference between 1\% and 0.1\% value at risk
#' @param qDiff 90\% quantile for mDiff
#' @param y Transformed data
#' @export

createExpertPrior <- function(m,q,mDiff,qDiff,y) {
	#Bayesian Analysis of Extreme Events with Threshold Estimation
	#Cibele N. Behrens, Hedibert F. Lopes and Dani Gamerman

	elicitedParams <- paramsFromExpert(m,q,mDiff,qDiff)

	a1 <- elicitedParams[1]
	b1 <- elicitedParams[2]
	a2 <- elicitedParams[3]
	b2 <- elicitedParams[4]
	p <- 0.95
	u <- thresholdFromPercentage(p,y)
	exceedances <- length(y[y>u])
	N <- length(y)
	p1 <- 0.01*N/exceedances
	p2 <- 0.001*N/exceedances

	expertPrior <- function(gp.parm) {
		sigma <- gp.parm[2]
		ksi <- gp.parm[3]
		sigmaksi.prior <- log(u+sigma/ksi*(p1^-ksi-1))*(a1-1) +
		  -b1*(u+sigma/ksi*(p1^-ksi-1)) +
		  log(sigma/ksi*(p2^-ksi - p1^-ksi))*(a2-1) +
		  -b2*(sigma/ksi*(p2^-ksi - p1^-ksi)) +
		  log(abs(-sigma/ksi^2 *
		            ((p1*p2)^-ksi*(log(p2)-log(p1))-p2^-ksi*log(p2)+p1^-ksi*log(p1))))
		  sigmaksi.prior
	}
	compiler::cmpfun(expertPrior)
}

#' Reference prior
#' @export

referencePrior <- function(gp.parm) {
	sigma <- gp.parm[2]
	ksi <- gp.parm[3]
	#Noninformative priors for the scale parameter in the generalized Pareto distribution
	#Sang Gil Kang
	sigma.prior <- log(1/(sigma*sqrt(1+ksi)*sqrt(1+2*ksi)))

	#Noninformative priors for the shape parameter in the generalized Pareto distribution
	#Sang Gil Kang, Dal Ho Kim, Woo Dong Lee
    ksi.prior <- log(1/(sigma*(1+ksi)*sqrt(1+2*ksi)))
    sigma.prior+ksi.prior
}

paramsFromExpert <- function(m,q,mDiff,qDiff) {

    outerFun <- function(median, q) {

	    innerFun <- function(shape) {
	        scale <- median/shape * (3*shape+0.2)/(3*shape-0.8)

	        pgamma(q, shape,,scale) - 0.90
	    }

	    tmp <- uniroot(innerFun, lower=1, upper=100000)

	    myshape <- tmp$root
	    myscale <- median/myshape * (3*myshape+0.2)/(3*myshape-0.8)

	    c(myshape, 1/myscale) #returning rate instead of scale
    }

    firstParams <- outerFun(m, q)
    a1 <- firstParams[1]
    b1 <- firstParams[2]

    secondParams <- outerFun(mDiff, qDiff)
    a2 <- secondParams[1]
    b2 <- secondParams[2]

    c(a1,b1,a2,b2)
}