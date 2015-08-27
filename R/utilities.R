#' Utility function for transforming stock data
#' @param orgData Non-transformed data
#' @export

transformStockData <- function(orgData) {
    -1+orgData[2:NROW(orgData)]/orgData[1:NROW(orgData)-1]
}

histAndPlot <- function(distFunc, y, parm, breaks=200) {
	#distFunc takes one argument
	threshold <- parm[length(parm)-1]
	y <- sort(y)
	h <- hist(y,breaks=breaks,xlim=c(min(y),max(y)), xlab='Negative Return', main=NA)

	body.index <- y<=threshold
  	len <- sum(body.index)+1

	xfit <- y
	yfitOrg <- distFunc(xfit)
	yfit <- yfitOrg*diff(h$mids[1:2])*length(y)
	lines(xfit[1:len],yfit[1:len],type='l', lty=3, lwd=2)
	lines(xfit[len:length(xfit)],yfit[len:length(yfit)],type='l', lty=1, lwd=2)
}

distFromModel <- function(parm,model,data) {
	function(x) {
		r <- rep(0,length(x))
		for(i in 1:length(x)) {
			d <- data
			d$y <- x[i]
			M <- model(parm,d)
			r[i] <- exp(-M[['Dev']]/2)*data$Ntotal/data$N
		}
		r
	}
}

estDistAndPlot <- function(distFunc, y, parm) {
	threshold <- parm[length(parm)-1]
	body.index <- y<=threshold

	plot(y,estimateDist(y)(y),type='l')
	lines(y[body.index],distFunc(y[body.index]),col='blue',lwd=2)
	lines(y[!body.index],distFunc(y[!body.index]),col='red',lwd=2)
}

estimateDist <- function(y) {
	y <- sort(y)
	length.y <- length(y)
	#l <- quantile(y[2:length.y]-y[1:{length.y-1}], probs=0.99, names=FALSE)
	l <- mean(y[2:length.y]-y[1:{length.y-1}])*2
	area <- length.y*l*2
	d <- function(x) {
		r <- rep(0,length(x))
		for(i in 1:length(x))
			r[i] <- (sum(abs(y-x[i])<l)-1)/area
		r
	}
	d
}

#' Utility function for calculating percentile from threshold
#' @param u Threshold
#' @param y Transformed data
#' @export

percentageFromThreshold <- function(u,y) {
	k <- sum(y<=u)
	k/length(y)
}

#' Utility function for calculating threshold from percentile
#' @param p Percentile
#' @param y Transformed data
#' @export

thresholdFromPercentage <- function(p,y) {
	y <- sort(y)
	if(p==0) return(y[1])
	k <- ceiling(p*length(y))
	y[k]
}

valueAtRisk <- function(p, gp.parm, factor, timeHorizon=1) {
	threshold <- gp.parm[1]
	sigma <- gp.parm[2]
	ksi <- gp.parm[3]
	VaR <- threshold - sigma/ksi * (1-(1/factor*(1-p))^(-ksi))
	VaR*timeHorizon^ksi
}
#valueAtRisk <- cmpfun(valueAtRisk)

expectedShortfall <- function(VaR, gp.parm) {
	threshold <- gp.parm[1]
	sigma <- gp.parm[2]
	ksi <- gp.parm[3]
	if(ksi <= 0 || ksi >= 1) return(NA)
	(VaR+sigma-ksi*threshold)/(1-ksi)
}
#expectedShortfall <- cmpfun(expectedShortfall)

#Helpers for accessing rugarch GH distribution
dgh <- function(x,param,log=FALSE) {
	rugarch:::dgh(x,param[3],param[4],param[2],param[1],param[5],log=log)
}
rgh <- function(n,param) {
	rugarch:::rgh(n,param[3],param[4],param[2],param[1],param[5])
}
pgh <- function(q,param) {
	Integral <- integrate(dgh, q, Inf, param, stop.on.error = TRUE, rel.tol = .Machine$double.eps^0.5)
	1-Integral$value
}
qgh <- function(p,param) {
	rugarch:::qgh(p,param[3],param[4],param[2],param[1],param[5])
}

dgp <- function(x, gp.parm, log=FALSE) {
	#only x can be vector
	mu <- gp.parm[1]
	sigma <- gp.parm[2]
	ksi <-  gp.parm[3]
	if(sigma <= 0) stop("The sigma parameter must be positive.")
	N <- length(x)
	inside <- x >= mu & (ksi >= 0 | x <= mu - sigma/ksi)
	z <- (x[inside] - mu) / sigma
	dens <- rep(-Inf, N)
	if(ksi == 0) {
		dens[inside] <- -z - log(sigma)
	} else {
		dens[inside] <- log(1/sigma) + log(1 + ksi * z) * (-1/ksi - 1)
	}
	if(log == FALSE) dens <- exp(dens)
	return(dens)
}
#dgp <- cmpfun(dgp)

pgp <- function(x,gp.parm) {
	#only handles one x
	mu <- gp.parm[1]
	sigma <- gp.parm[2]
	ksi <- gp.parm[3]
	ub <- Inf
	if(ksi < 0) ub <- mu-sigma/ksi
	if(x < mu) return(0)
	if(x > ub) return(1)
	if(ksi == 0) {
		p <- 1-exp(-(x-mu)/sigma)
	}
	else {
		p <- 1-(1+ksi*(x-mu)/sigma)^(-1/ksi)
	}
	p
}
#pgp <- cmpfun(pgp)

qgp <- function(p,gp.parm) {
	as.vector(fExtremes::qgpd(p,gp.parm[3],gp.parm[1],gp.parm[2]))
}

#' Function for mean residual life plot
#' @param transformedData Transformed data
#' @export

MRLP <- function(transformedData) {
	transformedData <- sort(transformedData)
	threshold <- transformedData
	N <- length(threshold)

	meanExcess <- rep(0,N)
	pVector <- rep(0,N)

	for(i in 1:N) {
	exceedances <- transformedData[transformedData > threshold[i]]
	numberOfExceedances <- length(exceedances)
	meanExcess[i] <- 1/numberOfExceedances * sum(exceedances-threshold[i])
	pVector[i] <- (N-i)/N
	}
	pVector <<- pVector
	meanExcess <<- meanExcess

	labelV <- c(0.9, 0.99)*N

	plot(threshold, meanExcess, xaxt='n',
	   #xlim=c(-0.05, transformedData[length(transformedData)-1]),
	   #ylim=c(0,max(meanExcess[is.finite(meanExcess)])*0.7),
	   xlab="Threshold",
	   ylab="Mean Excess")
	axis(1, at=threshold[labelV],labels=c("90%", "99%"))
	abline(v=transformedData[length(transformedData)*0.95],lty=3, lwd=2) # vertical line  

}

#' Fixed-threshold GP fitted using MLE
#' @param y Transformed data
#' @param gp.p Percentile
#' @param gp.u Threshold (optional, percentile ignored if defined)
#' @export

fitGP_MLE <- function(y,gp.p,gp.u=-1) {

	y <- sort(y)

	if(gp.u == -1) gp.u <- y[length(y)*gp.p]

	y <- y[y>gp.u]

	dgp <- function(x) {
		gp.parm <- c(gp.u,x)
		LL <- try(sum(dgpd(y,gp.parm[1],gp.parm[2],gp.parm[3],log=TRUE)), silent=TRUE)
		if(inherits(LL, 'try-error')) return(Inf)
		-LL
	}

	iv <- c(0.01,0.1)

	res <- optim(
		iv,
		dgp,
		method='Nelder-Mead',
		control=list(trace=TRUE,maxit=10000))

	if(res$convergence != 0) stop('MLE didnt converge, ', res$message)

	gp.parm <- c(gp.u,res$par)

	histAndPlot(function(x) {
		dgpd(x,gp.parm[1],gp.parm[2],gp.parm[3],log=FALSE)
	},y,gp.parm[1:2], breaks=30)

	gp.parm
}