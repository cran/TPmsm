TPCheck <- function(object, s, t) {
	if ( !is.survTP(object) ) stop("Argument 'object' must be of class 'survTP'")
	if ( !is.numeric(s) ) return("Argument 's' is not numeric")
	if ( !is.numeric(t) ) return("Argument 't' is not numeric")
	if ( !(0 <= s & s <= t) ) return("'s' and 't' must be positive, and s <= t")
	return(NULL)
}

TPCCheck <- function(object, s, t, x) {
	Message <- TPCheck(object, s, t)
	if ( !is.null(Message) ) return(Message)
	if ( !is.numeric(x) ) return("Argument 'x' is not numeric")
	return(NULL)
}

TPWindowCheck <- function(h, nh, ncv, window) {
	if ( !is.numeric(h) ) return("Argument 'h' must be numeric")
	if (length(h) < 1 | length(h) > 2) return("Argument 'h' length must be between 1 and 2")
	if ( any(h <= 0) ) return("Argument 'h' must be greater than 0")
	if ( !( is.numeric(nh) | is.integer(nh) ) ) return("Argument 'nh' must be numeric or integer")
	if (nh <= 1) return("Argument 'nh' must be greater than 1")
	if ( !( is.numeric(ncv) | is.integer(ncv) ) ) return("Argument 'ncv' must be numeric or integer")
	if (ncv < 10) return("Argument 'ncv' must be greater or equal than 10")
	window0 <- c("normal", "epanech", "biweight", "triweight")
	window1 <- c(window0, "box")
	window2 <- c("asymmetric", "symmetric")
	window3 <- c(window1, "tricube", "triangular", "cosine", window2)
	if ( !( window %in% window3 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight', 'triweight', 'box', 'tricube', 'triangular', 'cosine', 'asymmetric' or 'symmetric'")
	return(NULL)
}

TPCWindowCheck <- function(bw, window, method.weights) {
	if ( !( is.character(bw) | is.numeric(bw) ) ) return("Argument 'bw' must be either a character string or a numeric vector")
	if ( is.character(bw) ) {
		if ( !exists(bw, mode="function") ) return( paste("could not find function '", bw, "'", sep="") )
	}
	window0 <- c("normal", "epanech", "biweight", "triweight")
	window1 <- c(window0, "box")
	window2 <- c("asymmetric", "symmetric")
	window3 <- c(window1, "tricube", "triangular", "cosine", window2)
	if (bw %in% c("ALbw", "CVbw", "PBbw") & !( window %in% window0 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight' or 'triweight'")
	else if (bw == "dpik" & !( window %in% window1 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight', 'triweight' or 'box'")
	else if (bw == "NNEbw" & !( window %in% window2 ) ) return("Argument 'window' must be one of 'asymmetric' or 'symmetric'")
	if ( !( window %in% window3 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight', 'triweight', 'box', 'tricube', 'triangular', 'cosine', 'asymmetric' or 'symmetric'")
	if ( !method.weights %in% c("NW", "LL") ) return("Argument 'weights' must be one of 'NW' or 'LL'")
	return(NULL)
}

StateCheck <- function(state.names) {
	if (length(state.names) != 3) return("Argument 'state.names' length must be equal to 3")
	if ( length(state.names) != length( unique(state.names) ) ) return("Argument 'state.names' must be unique")
	return(NULL)
}

BootCheck <- function(conf, n.boot, conf.level, method.boot) {
	if ( !is.logical(conf) ) return("Argument 'conf' must be logical")
	if ( !( is.numeric(n.boot) | is.integer(n.boot) ) ) return("Argument 'n.boot' must be numeric or integer")
	if (n.boot <= 1) return("Argument 'n.boot' must be greater than 1")
	if ( !is.numeric(conf.level) ) return("Argument 'conf.level' is not numeric")
	if (conf.level < 0 | conf.level > 1) return("Argument 'conf.level' must be between 0 and 1")
	if ( !( method.boot %in% c("percentile", "basic") ) ) return("Argument 'method.boot' must be one of 'percentile' or 'basic'")
	return(NULL)
}

CvalCheck <- function(boot.cv) {
	if ( !is.logical(boot.cv) ) return("Argument 'boot.cv' must be logical")
	return(NULL)
}

TPStateBootCheck <- function(object, s, t, state.names, conf, n.boot, conf.level, method.boot) {
	Message <- TPCheck(object, s, t)
	if ( !is.null(Message) ) return(Message)
	Message <- StateCheck(state.names)
	if ( !is.null(Message) ) return(Message)
	Message <- BootCheck(conf, n.boot, conf.level, method.boot)
	return(Message)
}

TPWindowStateBootCvalCheck <- function(object, s, t, h, nh, ncv, window, state.names, conf, n.boot, conf.level, method.boot, boot.cv) {
	Message <- TPCheck(object, s, t)
	if ( !is.null(Message) ) return(Message)
	Message <- TPWindowCheck(h, nh, ncv, window)
	if ( !is.null(Message) ) return(Message)
	Message <- StateCheck(state.names)
	if ( !is.null(Message) ) return(Message)
	Message <- BootCheck(conf, n.boot, conf.level, method.boot)
	if ( !is.null(Message) ) return(Message)
	Message <- CvalCheck(boot.cv)
	return(Message)
}

TPCWindowStateBootCheck <- function(object, s, t, x, bw, window, method.weights, state.names, conf, n.boot, conf.level, method.boot) {
	Message <- TPCCheck(object, s, t, x)
	if ( !is.null(Message) ) return(Message)
	Message <- TPCWindowCheck(bw, window, method.weights)
	if ( !is.null(Message) ) return(Message)
	Message <- StateCheck(state.names)
	if ( !is.null(Message) ) return(Message)
	Message <- BootCheck(conf, n.boot, conf.level, method.boot)
	return(Message)
}
