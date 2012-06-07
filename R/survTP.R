survTP <- function(time1, event1, Stime, event, covariate=NULL) {
	if ( missing(time1) ) stop("Argument 'time1' is missing, with no default")
	if ( missing(event1) ) stop("Argument 'event1' is missing, with no default")
	if ( missing(Stime) ) stop("Argument 'Stime' is missing, with no default")
	if ( missing(event) ) stop("Argument 'event' is missing, with no default")
	if ( !is.numeric(time1) ) stop("Argument 'time1' is not numeric")
	if ( !( is.logical(event1) | is.numeric(event1) ) ) stop("Argument event1 must be logical or numeric")
	if ( !is.numeric(Stime) ) stop("Argument 'Stime' is not numeric")
	if ( !( is.logical(event) | is.numeric(event) ) ) stop("Argument event must be logical or numeric")
	len <- length(time1)
	if ( is.null(covariate) ) {
		if ( len != length(event1) | len != length(Stime) | len != length(event) ) stop("Arguments 'time1', 'event1', 'Stime' and 'event' must have the same length")
	} else {
		if ( !is.numeric(covariate) ) stop("Argument 'covariate' is not numeric")
		if ( len != length(event1) | len != length(Stime) | len != length(event) | len != length(covariate) ) stop("Arguments 'time1', 'event1', 'Stime', 'event' and 'covariate' must have the same length")
	}
	if ( any( (event1 != 0 & event1 != 1) | (event1 != FALSE & event1 != TRUE) ) ) stop("Argument 'event1' must be 0 or 1 if numeric and TRUE or FALSE if logical")
	if ( any( (event != 0 & event != 1) | (event != FALSE & event != TRUE) ) ) stop("Argument 'event' must be 0 or 1 if numeric and TRUE or FALSE if logical")
	if ( any(time1 < 0 | Stime < 0) ) stop("Arguments 'time1' and 'Stime' must be greater than 0")
	if ( any(Stime < time1) ) stop("Argument 'Stime' must be greater or equal to argument 'time1'")
	if ( any(!event1 & Stime != time1) ) stop("Arguments 'Stime' and 'time1' must be equal when argument 'event1' equals 0 or FALSE")
	if ( any(!event1 & event) ) stop("Argument 'event' must be equal to 0 or FALSE when argument 'event1' equals 0 or FALSE")
	if ( any(time1 == Stime & event1 & !event) ) stop("When arguments 'Stime' and 'time1' are equal and argument 'event1' equals 1 or TRUE, argument 'event' must equal 1 or TRUE")
	object <- vector(mode="list", length=1)
	if ( is.null(covariate) ) object[[1]] <- data.frame( "time1"=as.double(time1), "event1"=as.integer(event1), "Stime"=as.double(Stime), "event"=as.integer(event) )
	else object[[1]] <- data.frame( "time1"=as.double(time1), "event1"=as.integer(event1), "Stime"=as.double(Stime), "event"=as.integer(event), "covariate"=as.double(covariate) )
	class(object) <- "survTP"
	return(object)
}

is.survTP <- function(object) inherits(object, "survTP")
