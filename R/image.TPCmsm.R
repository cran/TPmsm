image.TPCmsm <- function(x, plot.type="tc", tr.choice, xlim, ylim, zlim, col, xlab, ylab, main, sub, key.title, key.axes, las=1, conf.int=FALSE, legend=TRUE, curvlab, contour=TRUE, nlevels, levels, ...) {
	if ( !inherits(x, "TPCmsm") ) stop("'x' must be of class 'TPCmsm'")
	if ( !( plot.type %in% c("tc", "ct") ) ) stop("Argument 'plot.type' must be one of 'tc' or 'ct'")
	if ( missing(tr.choice) ) tr.choice <- dimnames(x$est)[[3]]
	lt <- length(tr.choice)
	if (sum( tr.choice %in% dimnames(x$est)[[3]] ) != lt) stop("Argument 'tr.choice' and possible transitions must match")
	if ( anyDuplicated(tr.choice) ) stop("Argument 'tr.choice' must be unique")
	require(grDevices)
	require(graphics)
	itr <- match( tr.choice, dimnames(x$est)[[3]] )
	if ( missing(main) ) main <- ""
	if ( missing(sub) ) sub <- ""
	if ( missing(zlim) ) zlim <- c(0, 1)
	if ( missing(nlevels) ) nlevels <- 20
	if ( missing(levels) ) levels <- pretty(zlim, nlevels)
	if ( missing(col) ) col <- heat.colors(nlevels)[nlevels:1]
	if ( missing(curvlab) ) curvlab <- tr.choice
	mat <- matrix(nrow=2*conf.int+1, ncol=lt+1)
	par.orig <- par( c("mar", "las", "mfrow", "new", "cex") )
	on.exit( par(par.orig) )
	if (conf.int) {
		mat[1,1:lt] <- (lt*2+1):(lt*3)
		mat[2,1:lt] <- 1:lt
		mat[3,1:lt] <- (lt+1):(lt*2)
	} else mat[1,1:lt] <- 1:lt
	mat[,ncol(mat)] <- rep( (2*conf.int+1)*lt+1, nrow(mat) )
	mar.orig <- par.orig$mar
	if (length(tr.choice) < 2 & !conf.int) w <- (3 + mar.orig[2L]) * par("csi") * 2.54
	else w <- (3 + mar.orig[2L]) * par("csi") * 2
	layout( mat, widths=c( rep(1, lt), lcm(w) ) )
	mar <- mar.orig
	mar[4L] <- 1
	par(las=las, mar=mar)
	if (plot.type == "tc") {
		if ( missing(xlab) ) xlab <- "Time"
		if ( missing(ylab) ) ylab <- "Covariate"
		if ( missing(xlim) ) xlim <- c(x$time[1], x$time[length(x$time)])
		if ( missing(ylim) ) ylim <- c(x$covariate[1], x$covariate[length(x$covariate)])
		for ( i in seq_len(lt) ) {
			image(x=x$time, y=x$covariate, z=x$est[,,itr[i]], xlim=xlim, ylim=ylim, zlim=zlim, col=col, main="", sub="", xlab="", ylab="", breaks=levels, ...)
			if (contour) contour(x=x$time, y=x$covariate, z=x$est[,,itr[i]], nlevels=nlevels, levels=levels, xlim=xlim, ylim=ylim, zlim=zlim, axes=FALSE, col=grey(0.4), add=TRUE)
			if (legend) title(main=curvlab[i], sub="", xlab="", ylab="", ...)
		}
	} else if (plot.type == "ct") {
		if ( missing(xlab) ) xlab <- "Covariate"
		if ( missing(ylab) ) ylab <- "Time"
		if ( missing(xlim) ) xlim <- c(x$covariate[1], x$covariate[length(x$covariate)])
		if ( missing(ylim) ) ylim <- c(x$time[1], x$time[length(x$time)])
		for ( i in seq_len(lt) ) {
			image(x=x$covariate, y=x$time, z=t(x$est[,,itr[i]]), xlim=xlim, ylim=ylim, zlim=zlim, col=col, main="", sub="", xlab="", ylab="", breaks=levels, ...)
			if (contour) contour(x=x$covariate, y=x$time, z=t(x$est[,,itr[i]]), nlevels=nlevels, levels=levels, xlim=xlim, ylim=ylim, zlim=zlim, axes=FALSE, col=grey(0.4), add=TRUE)
			if (legend) title(main=curvlab[i], sub="", xlab="", ylab="", ...)
		}
	}
	if ( conf.int & !is.null(x$inf) & !is.null(x$sup) ) {
		if (plot.type == "tc") {
			for ( i in seq_len(lt) ) {
				image(x=x$time, y=x$covariate, z=x$inf[,,itr[i]], xlim=xlim, ylim=ylim, zlim=zlim, col=col, main="", sub="", xlab="", ylab="", breaks=levels, ...)
				if (contour) contour(x=x$time, y=x$covariate, z=x$inf[,,itr[i]], nlevels=nlevels, levels=levels, xlim=xlim, ylim=ylim, zlim=zlim, axes=FALSE, col=grey(0.4), add=TRUE)
			}
			for ( i in seq_len(lt) ) {
				image(x=x$time, y=x$covariate, z=x$sup[,,itr[i]], xlim=xlim, ylim=ylim, zlim=zlim, col=col, main="", sub="", xlab="", ylab="", breaks=levels, ...)
				if (contour) contour(x=x$time, y=x$covariate, z=x$sup[,,itr[i]], nlevels=nlevels, levels=levels, xlim=xlim, ylim=ylim, zlim=zlim, axes=FALSE, col=grey(0.4), add=TRUE)
			}
		} else if (plot.type == "ct") {
			for ( i in seq_len(lt) ) {
				image(x=x$covariate, y=x$time, z=t(x$inf[,,itr[i]]), xlim=xlim, ylim=ylim, zlim=zlim, col=col, main="", sub="", xlab="", ylab="", breaks=levels, ...)
				if (contour) contour(x=x$covariate, y=x$time, z=t(x$inf[,,itr[i]]), nlevels=nlevels, levels=levels, xlim=xlim, ylim=ylim, zlim=zlim, axes=FALSE, col=grey(0.4), add=TRUE)
			}
			for ( i in seq_len(lt) ) {
				image(x=x$covariate, y=x$time, z=t(x$sup[,,itr[i]]), xlim=xlim, ylim=ylim, zlim=zlim, col=col, main="", sub="", xlab="", ylab="", breaks=levels, ...)
				if (contour) contour(x=x$covariate, y=x$time, z=t(x$sup[,,itr[i]]), nlevels=nlevels, levels=levels, xlim=xlim, ylim=ylim, zlim=zlim, axes=FALSE, col=grey(0.4), add=TRUE)
			}
		}
	}
	mar <- par("mar")
	mar[4L] <- mar[2L]
	mar[2L] <- 1
	par(mar=mar)
	plot.new()
	plot.window(xlim=c(0, 1), ylim=range(levels), xaxs="i", yaxs="i")
	rect(0, levels[-length(levels)], 1, levels[-1L], col=col)
	if ( missing(key.axes) ) axis(4)
	else key.axes
	box()
	if ( missing(key.title) ) {
		if (length(tr.choice) < 2 & !conf.int) cex.main <- 0.8
		else cex.main <- 1
		title(main="Transition\nprobability", cex.main=cex.main)
	} else key.title
	par(par.orig)
	par(mfrow=c(1, 1), new=TRUE, cex=1.2)
	plot(xlim, ylim, main=main, sub=sub, xlab=xlab, ylab=ylab, axes=FALSE, type="n", bty="n", ...)
	invisible()
}
