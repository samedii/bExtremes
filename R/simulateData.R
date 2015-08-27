#' Creates expert prior from ellicitation
#' @param n Sample size, defaults to 1000
#' @export

simulateGHGP <- function(n=1000) {
	gh.parm <- c(0.002200946, 0.031815232, 15.241766213, -12.325859594, -3.336424312)
	gp.p <- 0.95
	gp.u <- qgh(gp.p,gh.parm)
	y <- rgh(n,gh.parm)
	index <- y>gp.u
	gp.sigma <- (1-gp.p)/dgh(gp.u,gh.parm,log=FALSE)
	gp.parm <- c(gp.u,gp.sigma,0.3)
	y[index] <- LaplacesDemon:::rgpd(sum(index),gp.parm[1],gp.parm[2],gp.parm[3])

	cat('Analytical values:\n')
	VaR <- valueAtRisk(0.99, gp.parm, 1-gp.p)
	VaR01 <- valueAtRisk(0.999, gp.parm, 1-gp.p)
	ES <- expectedShortfall(VaR, gp.parm)
	ES01 <- expectedShortfall(VaR01, gp.parm)

	cat('VaR0.01: ', VaR, '\n')
	cat('ES0.01: ', ES, '\n')
	cat('VaR0.001: ', VaR01, '\n')
	cat('ES0.001: ', ES01, '\n')

	y
}