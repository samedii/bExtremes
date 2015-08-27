#' Create fixed-threshold GP data object
#' @param model Model object
#' @param d Data object
#' @param ivAndCov Starting values and covariance
#' @export

fitMCMC <- function(model,d,ivAndCov) {

	iv <- LaplacesDemon::as.initial.values(ivAndCov)
	histAndPlot(distFromModel(iv,model,d),d$y,iv)

	L <- length(iv)

	if(d$dist == 'GH') {
		w <- c(0.000552087233911274, 0.00118877162812569, 0.00347019556493294,
			0.00948499325675634, 0.000487120914819678, 1.04061614433659e-05,
			3.33928564919788e-07)
	} else if(d$dist == 'GP') {
		w <- c(1,1,1,1)
	} else stop('Unrecognized distribution')

	iter = 1000
	Fit.adapted <- LaplacesDemon::LaplacesDemon(model, Data=d, Initial.Values=iv,
	    Covar=ivAndCov$Covar,
	    Iterations=iter, Status=50, Thinning=1,
	    Algorithm='AFSS', Specs=list(
	    	A=Inf,
	        B=NULL,
	        m=50, #interval at which heuristics are used
	        n=0,
	        w=w
	    ),
		Debug=list(DB.chol=TRUE, DB.eigen=TRUE, DB.MCSE=TRUE, DB.Model=FALSE))


	if(!LaplacesDemon::is.stationary(Fit.adapted)) {
		cat('adapted NOT STATIONARY, retrying...\n')

		iv.adapted <- LaplacesDemon::as.initial.values(Fit.adapted)
		try(histAndPlot(distFromModel(iv.adapted,model,d),d$y,iv.adapted))

		iter2 <- 5000
		Fit.adapted <- LaplacesDemon::LaplacesDemon(model, Data=d, iv.adapted,
		     Covar=Fit.adapted$Covar, Iterations=iter2, Status=50, Thinning=1,
		     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=Fit.adapted$Specs$m,
		     n=iter, w=Fit.adapted$CovarDHis[nrow(Fit.adapted$CovarDHis),]),
		     Debug=list(DB.chol=TRUE, DB.eigen=TRUE, DB.MCSE=TRUE, DB.Model=FALSE))
		iter <- iter+iter2

		if(!LaplacesDemon::is.stationary(Fit.adapted)) {
			iv.adapted <- LaplacesDemon::as.initial.values(Fit.adapted)
			histAndPlot(distFromModel(iv.adapted,model,d),d$y,iv.adapted)
			stop('adapted STILL NOT STATIONARY\n')
		}
	}

	iv.adapted <- LaplacesDemon::as.initial.values(Fit.adapted)
	mu.adapted <- Fit.adapted$Summary2[1:L,1]
	histAndPlot(distFromModel(mu.adapted,model,d),d$y,mu.adapted)

	Fit.short <- Fit.adapted
	attempts <- 0
	repeat{
		attempts <- attempts + 1
		iv.rep <- LaplacesDemon::as.initial.values(Fit.short)
		mu.rep <- Fit.short$Summary1[1:L,1]

		Fit.short.temp <- LaplacesDemon::LaplacesDemon(model, Data=d, iv.rep,
		     Covar=Fit.short$Covar*1.1, Iterations=10000, Status=500, Thinning=1,
		     Algorithm='IM', Specs=list(mu=mu.rep),
		     Debug=list(DB.chol=TRUE, DB.eigen=TRUE, DB.MCSE=TRUE, DB.Model=FALSE))

		cat('Acceptance rate: ', Fit.short.temp$Acceptance.Rate, '\n')
		if(Fit.short.temp$Acceptance.Rate > 0.15) {
			if(Fit.short.temp$Acceptance.Rate > 0.40) {
				Fit.short$Covar <- Fit.short$Covar*2
			}
			else {
				break
			}
		}
		else {
			Fit.short <- Fit.short.temp
		}

		if(attempts >= 5) {
			stop('short acceptance rate not good enough after multiple attempts')
		}
	}



	iv.short <- LaplacesDemon::as.initial.values(Fit.short)
	mu.short <- Fit.short$Summary2[1:L,1]
	histAndPlot(distFromModel(mu.short,model,d),d$y,mu.short)

	Fit.sampling <- LaplacesDemon::LaplacesDemon(model, Data=d, iv.short,
	     Covar=Fit.short$Covar*1.1,
	     Iterations=50000, Status=500,
	     Thinning=min(30,Fit.short.temp$Rec.Thinning),
	     Algorithm='IM', Specs=list(mu=mu.short),
	     Debug=list(DB.chol=TRUE, DB.eigen=TRUE, DB.MCSE=TRUE, DB.Model=FALSE))

	if(!LaplacesDemon::is.stationary(Fit.sampling)) {
		stop('sampling NOT STATIONARY\n')
	}

	if(!LaplacesDemon::is.appeased(Fit.sampling)) {
		Consort(Fit.sampling)
		stop('sampling NOT APPEASED\n')
	}

	mu.sampling <- Fit.sampling$Summary2[1:L,1]
	histAndPlot(distFromModel(mu.sampling,model,d),d$y,mu.sampling)

	Fit.sampling
}