# bExtremes
R package for extreme risk estimation using extreme value theory and Bayesian inference.

Automatic threshold weighting is done with body-tail models. Implemented body distributions include "GH" and "GP". The latter is recommended.

A fixed-threshold GP model is also included for comparisons (see functions createGPData and createGPModel).

## Usage

    loadNamespace('devtools')
    devtools::dev_mode(on=TRUE)
    
    #Requires LaplacesDemon package from github repo: samedii/LaplacesDemon
    attempt <- try(loadNamespace('LaplacesDemon'))
    if(inherits(attempt, 'try-error')) {
        devtools::install_github('samedii/LaplacesDemon')
    }
    
    attempt <- try(library(bExtremes))
    if(inherits(attempt, 'try-error')) {
        devtools::install_github('samedii/bExtremes')
    }
    library(bExtremes)
    
    #cat('Getting data...\n')
    #orgData <- fetchExampleStockData('BAC') #function unavailable, replace with your own
    
    #Transforming and flipping data (because right tail is modelled)
    #transformedData <- -transformStockData(orgData)
    
    #cat('Filtering volatility...\n')
    #transformedData <- filterVolatility(transformedData)
    
    cat('Simulating data...\n')
    transformedData <- simulateGHGP()
    
    cat('Creating prior...\n')
    #Expert prior (Guessing from historical data)
    sortedData <- sort(transformedData)
    ##Median for VaR at 0.01
    m <- sortedData[ceiling(length(sortedData)*(1-0.01))]
    ##90% quantile for VaR 0.01
    q <- m*1.2
    
    ##Median for difference between VaR at 0.001 and 0.01
    mDiff <- sortedData[ceiling(length(sortedData)*(1-0.001))]-m
    ##90% quantile of above
    qDiff <- mDiff*1.2
    
    prior <- createExpertPrior(m,q,mDiff,qDiff,transformedData)
    
    #Uninformed prior (reference prior)
    prior <- referencePrior
    
    #Calculate VaR and ES for confidence level p (preferably 0.99 or greater):
    p <- 0.99
    
    cat('Setting up model...\n')
    data <- createData(transformedData,p,dist='GP')
    model <- createModel(data,prior)
    
    cat('Generating IV...\n')
    ivAndCov <- generateIV(model,data,n=500) #It can be helpful to increase n
    
    cat('MCMC process...\n')
    fit <- fitMCMC(model,data,ivAndCov)
    
    #Information about final MCMC run and results
    LaplacesDemon::Consort(fit)
    
    #Summary of posterior samples
    fit$Summary2
    
    #Samples
    head(fit$Posterior2) #parameters
    head(fit$Monitor) #monitored variables (including VaR and ES)
    
    #Example plots
    hist(fit$Monitor[,7],breaks=100)
    LaplacesDemon::joint.density.plot(fit$Monitor[,7],fit$Monitor[,8])
    LaplacesDemon::joint.pr.plot(fit$Posterior2[,3],fit$Monitor[,7])
