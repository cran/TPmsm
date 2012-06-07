TransWidth <- function(object, bw, window, ...) {
	if ( is.numeric(bw) ) h <- bw
	else {
		if (bw =="dpik") {
			require(KernSmooth)
			h <- match.fun(bw)(x=object[[1]]$covariate, kernel=window, ...)
		}
		else if ( bw %in% c("ALbw", "CVbw", "PBbw") ) {
			type_kernel <- switch( window, "normal"="n", "epanech"="e", "biweight"="b", "triweight"="t" )
			if ( bw %in% c("ALbw", "PBbw") ) h <- match.fun(bw)(type_kernel=type_kernel, vec_data=object[[1]]$covariate, ...)
			else h <- match.fun(bw)(type_kernel=type_kernel, vec_data=object[[1]]$covariate, ...)$bw
		}
		else h <- match.fun(bw)(object[[1]]$covariate, ...)
	}
	return(h)
}
