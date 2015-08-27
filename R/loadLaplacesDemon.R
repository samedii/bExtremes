devtools::dev_mode(on=TRUE)
attempt <- try(loadNamespace('LaplacesDemon'))
if(inherits(attempt, 'try-error')) {
	install_github('samedii/LaplacesDemon')
}