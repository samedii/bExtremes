#' Volatility filter of transformed data
#' @param transformedData Transformed data
#' @export

filterVolatility <- function(transformedData) {
	spec = ugarchspec(
		variance.model = list(
			model = 'gjrGARCH',
			garchOrder = c(1,1)
			),
		mean.model = list(
			armaOrder = c(0,0),
			include.mean = FALSE
			),
		distribution.model = 'ghyp',
		start.pars = list(skew=0.2,shape=0.7,ghlambda=2)
		)
	fit = ugarchfit(
		data = transformedData,
		spec = spec,
		solver = 'solnp'
		)

	forc = ugarchforecast(fit, n.ahead=1)

	filteredData <- (as.vector(residuals(fit,standardize=TRUE))+fitted(forc))*sigma(forc)

	filteredData
}