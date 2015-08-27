generateGPIV <- function(n,d) {
	e <- .Machine$double.eps^0.5

	best.gp.parm <- rep(0,5)
	best.gp.LL <- -Inf

	for(i in 1:n) {

		sigma <- LaplacesDemon::interval(rnorm(1,0.1,0.02),e,Inf)
		ksi <- rnorm(1,0.01,0.1)

		body.gp.parm <- c(d$body.gp.u,sigma,ksi)

		try({
			LL <- sum(dgp(d$y,body.gp.parm,log=TRUE))
			if(LL > best.gp.LL) {
				best.gp.LL <- LL
				best.gp.parm <- body.gp.parm
			}
		})

	}

	if(!is.finite(best.gp.LL)) stop('Couldnt find a proper iv for GH')

	best.gp.parm[2:3]
}

generateGHIV <- function(n,d) {
	e <- .Machine$double.eps^0.5

	best.gh.parm <- rep(0,5)
	best.gh.LL <- -Inf

	for(i in 1:n) {
		gh.mu <- rnorm(1, 0, 0.01)
		gh.delta <- rgamma(1,10,500)#rnorm(1, 0.02, 0.02)
	    gh.alpha <- LaplacesDemon::rhalfnorm(1,0.1)#rnorm(1, 10, 3)
	    gh.beta <- runif(1,-abs(gh.alpha+e), abs(gh.alpha-e))#rnorm(1, 0, 5)
	    if(is.na(gh.beta)) gh.beta <- 0
	    gh.lambda <- rnorm(1, -2, 0.3)

		#gh.mu <- interval(gh.mu, -Inf, Inf)
		gh.delta <- LaplacesDemon::interval(gh.delta, e, Inf)
		gh.alpha <- LaplacesDemon::interval(gh.alpha, e, Inf)
		#gh.beta <- interval(gh.beta, -gh.alpha+1e-100, gh.alpha-1e-100)
		#gh.lambda <- interval(gh.lambda, -Inf, Inf)

		gh.parm <- c(gh.mu,gh.delta,gh.alpha,gh.beta,gh.lambda)

		try({
			LL <- sum(dgh(d$y,gh.parm,log=TRUE))
			if(LL > best.gh.LL) {
				best.gh.LL <- LL
				best.gh.parm <- gh.parm
			}
		})

	}

	if(!is.finite(best.gh.LL)) stop('Couldnt find a proper iv for GH')

	best.gh.parm
}

generateTailIV <- function(n,d,distFunc) {
	e <- .Machine$double.eps^0.5

	best.tail.parm <- c(0.05,0,0)
	best.LL <- -Inf

	improvements <- 0
	for(i in 1:n) {
		gp.u.center <- d$gp.u.ub-(d$gp.u.ub-d$gp.u.lb)*0.1
		gp.u <- LaplacesDemon::interval(rnorm(1,gp.u.center,0.1),d$gp.u.lb,d$gp.u.ub)
		#gp.sigma <- interval(rnorm(1,0.01,0.02),e,Inf)
		gp.ksi <- LaplacesDemon::interval(rnorm(1,0.5,1),e,1)

		tail.parm <- c(gp.u,gp.ksi)

		LL <- try(sum(distFunc(tail.parm)), silent=TRUE)

		if(!is.finite(LL) | inherits(LL, "try-error")) next

		if(LL > best.LL) {
			best.LL <- LL
			best.tail.parm <- tail.parm
			improvements <- improvements + 1
		}
	}

	if(!is.finite(best.LL)) stop('Couldnt find a proper iv for tail')

	best.tail.parm
}

#' Generate initial values and covariance matrix
#' @param model Model object
#' @param d Data object
#' @param n Number of guesses, defaults to 500
#' @export

generateIV <- function(model,d,n=500) {

	failFreeModel <- function(parm,d) {
		M <- try(model(parm,d),silent=TRUE)
		if(inherits(M, 'try-error')) return(list(LP=-Inf,LL=Inf,parm=parm))
		M
	}

	attempts <- 0
	repeat {
		attempts <- attempts +  1
		#Guess body.iv
		if(d$dist == 'GH') body.iv <- generateGHIV(n,d)
		else if(d$dist == 'GP') body.iv <- generateGPIV(n,d)
		else if(d$dist == 'NIG') stop('TODO')
		else stop('Unrecognized distribution')
		#Guess tail.iv
		tail.iv <- generateTailIV(n,d,function(tail.parm) {
				M <- try(model(c(body.iv,tail.parm),d),silent=TRUE)
				if(inherits(M, 'try-error')) return(-Inf)
				M$LP
			})
		histAndPlot(distFromModel(c(body.iv,tail.iv),model,d),d$y,c(body.iv,tail.iv))
		#Improve IV and estimate covariance
		Fit.la <- LaplacesDemon::LaplaceApproximation(Model=failFreeModel, c(body.iv,tail.iv), Data=d,
			Iterations=5000, Method='NM',
			Samples=d$N, CovEst='Hessian', sir=FALSE)

		histAndPlot(distFromModel(LaplacesDemon::as.initial.values(Fit.la),
			model,d),d$y,
			LaplacesDemon::as.initial.values(Fit.la))

		if(all(diag(Fit.la$Covar) == 1)) {
			if(attempts < 2) {
				cat('Failed to estimate covariance, retrying...\n')
			}
			else {
				cat('Failed to estimate covariance multiple times, trying guess and IM')

				if(d$dist == 'GH') {
					Fit.la$Covar <- structure(c(8.68953860504045e-07, 3.04538740278266e-07, -0.000641373326739513, 
						-0.00352409907933533, -3.78567958082721e-05, -6.46473671817658e-11, 
						3.48945900321482e-05, 3.04538740278266e-07, 2.94730664579168e-06, 
						6.77636876236097e-05, -0.00140073347066162, -0.000307215232859961, 
						1.6253865570465e-10, 3.67877360486246e-05, -0.000641373326739513, 
						6.77636876236097e-05, 2.22044604925031e-16, 2.9931703046162, 
						-0.0588204170968006, 3.16973066465142e-09, -0.0310548903855058, 
						-0.00352409907933533, -0.00140073347066162, 2.9931703046162, 
						16.4618184979521, 0.173794421986846, 3.23274062305456e-07, -0.160751522943159, 
						-3.78567958082721e-05, -0.000307215232859961, -0.0588204170968006, 
						0.173794421986846, 0.0360090035090007, -6.90054195781069e-09, 
						-0.00502994826475813, -6.46473671817658e-11, 1.6253865570465e-10, 
						3.16973066465142e-09, 3.23274062305456e-07, -6.90054195781069e-09, 
						3.4684493007183e-08, 2.04595954452649e-07, 3.48945900321482e-05, 
						3.67877360486246e-05, -0.0310548903855058, -0.160751522943159, 
						-0.00502994826475813, 2.04595954452649e-07, 0.0108146451537264
						), .Dim = c(7L, 7L), .Dimnames = list(c("gh.mu", "gh.delta", 
						"gh.alpha", "gh.beta", "gh.lambda", "gp.u", "gp.ksi"), c("gh.mu", 
						"gh.delta", "gh.alpha", "gh.beta", "gh.lambda", "gp.u", "gp.ksi"
						)))
				} else if(d$dist == 'GP') {
					Fit.la$Covar <- structure(c(3.32603710300194e-10, 3.32778916941703e-10, 3.32884741700801e-10, 
						3.32664361413157e-10, 3.32778916941703e-10, 3.32958253965702e-10, 
						3.33061611568457e-10, 3.32841439217669e-10, 3.32884741700801e-10, 
						3.33061611568457e-10, 3.33168549923906e-10, 3.32946628690475e-10, 
						3.32664361413157e-10, 3.32841439217669e-10, 3.32946628690475e-10, 
						3.32729476675301e-10), .Dim = c(4L, 4L), .Dimnames = list(c("body.gp.sigma", 
						"body.gp.ksi", "gp.u", "gp.ksi"), c("body.gp.sigma", "body.gp.ksi", 
						"gp.u", "gp.ksi")))
				} else stop('Unrecognized distribution')

				iv.la <- LaplacesDemon::as.initial.values(Fit.la)
				Fit.covar <- LaplacesDemon::LaplacesDemon(model, Data=d, iv.la,
	     			Covar=Fit.la$Covar*1.1, Iterations=5000, Status=500, Thinning=1,
	     			Algorithm='IM', Specs=list(mu=iv.la),
	     			Debug=list(DB.chol=TRUE, DB.eigen=TRUE, DB.MCSE=TRUE, DB.Model=FALSE))

				Fit <- Fit.la
				Fit$Covar <- Fit.covar$Covar
				return(Fit)
			}
		}
		else {
			Fit.la$Covar <- Fit.la$Covar*2.38^2/d$parm.length
			return(Fit.la)
		}
	}

}