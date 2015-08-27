#' Create body-tail model
#' @param data Data object from "createData"
#' @param gp.prior Tail-GP prior
#' @export

createModel <- function(data, gp.prior) {

	if(data$dist == 'GH') {
		dbody <- function(x,gh.parm,log=FALSE) {
			dgh(x,gh.parm,log)
		}
		pbody <- function(x,body.parm) {
			Integral <- integrate(dbody,
			x,Inf,
			body.parm,
			log=FALSE,
			stop.on.error=TRUE,rel.tol=.Machine$double.eps^0.5)
			1-Integral$value
		}
		body.prior <- function(body.parm) {
			#Maybe these are overly informed?
		    prior.mu <- dnorm(body.parm[1], 0, 1e-02, log=TRUE)
			prior.delta <- dgamma(body.parm[2], 10, 500, log=TRUE)
	    	prior.alpha <- LaplacesDemon::dhalfnorm(body.parm[3], 0.1, log=TRUE)
	    	prior.beta <- dunif(body.parm[4], -abs(body.parm[3]), abs(body.parm[3]), log=TRUE)
	    	prior.lambda <- dnorm(body.parm[5], -2, 3, log=TRUE)
			prior.mu + prior.delta + prior.alpha + prior.beta + prior.lambda
		}
	} else if(data$dist == 'GP') {
		dbody <- function(x,body.parm,log=FALSE) {
			if(body.parm[2] <= -1) stop('body.gp.ksi <= -1')
			db <- dgp(x,c(data$body.gp.u,body.parm[1],body.parm[2]),log)
			if(log) {
				db <- db + base::log(1-data$body.gp.p)
			} else {
				db <- db*(1-data$body.gp.p)
			}
			db
		}
		pbody <- function(x,body.parm) {
			pgp(x,c(data$body.gp.u,body.parm))*(1-data$body.gp.p)+data$body.gp.p
		}
		body.prior <- function(body.parm) {
			prior.sigma <- dnorm(body.parm[1],0,10)
			prior.ksi <- dnorm(body.parm[2],0,10)
			prior.sigma + prior.ksi
		}
	} else if(data$dist == 'NIG') {
		stop('TODO')
	} else {
		stop('Unrecognized distribution')
	}
	dbody <- compiler::cmpfun(dbody)
	pbody <- compiler::cmpfun(pbody)
	body.prior <- compiler::cmpfun(body.prior)

	gp.lb <- c(data$gp.u.lb,0,0)
	gp.ub <- c(data$gp.u.ub,Inf,1)

	dtail <- function(x,gp.parm,log=FALSE) {
		dgp(x,gp.parm,log)
	}
	dtail <- compiler::cmpfun(dtail)

	body.index <- 1:(data$parm.length-2)
	tail.index <- (data$parm.length-1):data$parm.length

	e <- .Machine$double.eps^0.5

	model <- function(parm, d) {

		body.parm <- parm[body.index]
		tail.parm <- parm[tail.index]

		#<body>-gptail distribution
		LL <- 0

	    ##Body
	    body.dataset.index <- d$y<=tail.parm[1]
	    body.dataset <- d$y[body.dataset.index]
	    if(length(body.dataset) > 0) {
			LL <- sum(dbody(body.dataset,body.parm,log=TRUE))

	        if(!is.finite(LL)) {
	            stop('Infinite in body. parm: ', parm, '\n')
	        }
	    }

	    ##Generalized Pareto
	    ###gp.p.cdf (used to adjust height of dist)
		gp.p.cdf <- pbody(tail.parm[1],body.parm)
    	gp.p.ex <- 1-sum(data$y>tail.parm[1])/d$Ntotal

    	if(abs(gp.p.cdf-gp.p.ex)>0.05) {
    		stop('Big gp.p diff, cdf: ', gp.p.cdf, ' vs ex:', gp.p.ex, '\n')
    	}

    	###gp.sigma
		gp.sigma <- (1-gp.p.cdf)/dbody(tail.parm[1],body.parm,log=FALSE)
		gp.parm <- c(tail.parm[1],gp.sigma,tail.parm[2])

	    if(any(gp.lb >= gp.parm | gp.parm >= gp.ub)) {
			stop('gp bounds not satisfied by: ', gp.parm)
		}

	    gp.dataset.index <- !body.dataset.index
	    gp.dataset <- d$y[gp.dataset.index]
	    if(length(gp.dataset) > 0) {
	    	LL.gp.one <- sum(dtail(gp.dataset,gp.parm,log=TRUE))

	    	LL <- LL + LL.gp.one + log(1-gp.p.cdf)*length(gp.dataset)

	        if(!is.finite(LL))
	            stop('Infinite in tail. param: ', parm, '\n')
	    }

		#Likelihood
		LL.prior <- body.prior(body.parm) + gp.prior(gp.parm)
	    LP <- LL + LL.prior

		#Value at Risk and Expected Shortfall
		factor <- 1-gp.p.cdf
		VaR <- valueAtRisk(d$p, gp.parm, factor)
		ES <- expectedShortfall(VaR, gp.parm)

	    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,LL.prior,gp.sigma,gp.p.cdf,gp.p.ex,gp.p.cdf-gp.p.ex,VaR,ES),
	                    yhat=NULL,
	                    parm=parm)
	    return(Modelout)
	}
	compiler::cmpfun(model)
}

#' Create fixed-threshold GP model
#' @param gp.prior GP prior
#' @export

createGPModel <- function(gp.prior) {

	gp.lb <- c(-Inf,0,0)
	gp.ub <- c(Inf,Inf,1)

	gpModel <- function(parm,d) {
		gp.parm <- c(d$gp.u,parm)

	    if(any(gp.lb >= gp.parm | gp.parm >= gp.ub)) {
			stop('gp bounds not satisfied')
		}

		LL <- sum(dgpd(d$y,gp.parm[1],gp.parm[2],gp.parm[3],log=TRUE))
		LL.prior <- gp.prior(gp.parm)
		LP <- LL + LL.prior

		#Value at Risk and Expected Shortfall
		factor <- 1-d$gp.p
		VaR <- valueAtRisk(d$p, gp.parm, factor)
		ES <- expectedShortfall(VaR, gp.parm)

		Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,LL.prior,VaR,ES),
					    yhat=NULL,#LaplacesDemon::rgpd(d$N,gp.parm[1],gp.parm[2],gp.parm[3]),
					    parm=parm)
		   return(Modelout)
	}
	compiler::cmpfun(gpModel)
}