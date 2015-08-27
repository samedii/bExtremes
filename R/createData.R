#' Create body-tail data object
#' @param y Transformed data
#' @param p Confidence level for risk measures
#' @param dist Body distribution for the model
#' @export

createData <- function(y, p, dist='GH') {

	y <- sort(y)

	Ntotal <- length(y)

	gp.p.lb <- 0.9
	if(p < 0.95) stop('p should not be smaller than 0.95')
	if(p < 0.99) warning('sanity of model becomes less sound with p<0.99')
	gp.p.ub <- p
	if(p > 0.99) gp.p.lb <- 0.99
	gp.u.lb <- thresholdFromPercentage(gp.p.lb,y)
	gp.u.ub <- thresholdFromPercentage(gp.p.ub,y)

	body.data <- list()
	if(dist == 'GH') {
		body.names <- list(
		    gh.mu=0,
		    gh.delta=0,
		    gh.alpha=0,
		    gh.beta=0,
		    gh.lambda=0
			)
	} else if(dist == 'GP') {

		body.gp.p <- 0.85
		body.gp.u <- thresholdFromPercentage(body.gp.p,y)
		y <- y[y>body.gp.u]
		body.data <- c(body.data,list(body.gp.p=body.gp.p,body.gp.u=body.gp.u))

		body.names <- list(
			body.gp.sigma=0,
			body.gp.ksi=0
			)
	} else if(dist == 'NIG') {
		stop('TODO')
	} else {
		stop('Unrecognized distribution')
	}

	mon.names <- list('LP','LL.prior','gp.sigma',
		'gp.p.cdf','gp.p.ex','gp.p.diff',
		'VaR','ES')

	tail.names <- list(
	    gp.u=0,
	    gp.ksi=0
	    )
	parm.names <- LaplacesDemon::as.parm.names(c(body.names,tail.names))

	r <- list(
		N=length(y),
		Ntotal=Ntotal,
		y=y,
		mon.names=mon.names,
		parm.names=parm.names,
		parm.length=length(parm.names),
		dist=dist,
		p=p,
		gp.p.lb=gp.p.lb, gp.p.ub=gp.p.ub,
		gp.u.lb=gp.u.lb, gp.u.ub=gp.u.ub
		)

	c(r,body.data)
}

#' Create fixed-threshold GP data object
#' @param y Transformed data
#' @param p Confidence level for risk measures
#' @param gp.p Percentile for threshold
#' @export

createGPData <- function(y,p,gp.p) {
	Ntotal <- length(y)

	e <- .Machine$double.eps^0.5

	y <- sort(y)
	gp.u <- thresholdFromPercentage(gp.p,y)-e

	y <- y[y>gp.u]

	mon.names <- list('LP','LL.prior','VaR','ES')

	parm.names <- LaplacesDemon::as.parm.names(list(
	    gp.sigma=0,
	    gp.ksi=0
	    ))

	r <- list(
		N=length(y),
		Ntotal=Ntotal,
		y=y,
		mon.names=mon.names,
		parm.names=parm.names,
		gp.u=gp.u,
		gp.p=gp.p,
		p=p,
		dist='fixed'
		)

	r
}