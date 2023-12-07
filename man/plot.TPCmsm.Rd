\encoding{UTF-8}
\name{plot.TPCmsm}
\alias{plot.TPCmsm}
\title{plot method for a TPCmsm object}

\description{
plot method for an object of class \sQuote{TPCmsm}.
It draws the estimated transition probabilities in a basic scatterplot.
}

\usage{\S3method{plot}{TPCmsm}(x, plot.type="t", tr.choice, xlab, ylab, col, lty, xlim, ylim,
conf.int=FALSE, ci.col, ci.lty, legend=TRUE, legend.pos, curvlab,
legend.bty="n", ...)}

\arguments{
	\item{x}{An object of class \sQuote{TPCmsm}.}
	\item{plot.type}{A character string specifying the type of plot.
		If \dQuote{t} the scatterplot of transition probability versus time is plotted.
		If \dQuote{c} the scatterplot of transition probability versus covariate is plotted.
	}
	\item{tr.choice}{Character vector of the form \sQuote{c(\dQuote{from to}, \dQuote{from to})}
		specifying which transitions should be plotted.
		Default, all the transition probabilities are plotted.
	}
	\item{xlab}{x-axis label.}
	\item{ylab}{y-axis label.}
	\item{col}{Vector of colour.}
	\item{lty}{Vector of line type. Default is 1:number of transitions.}
	\item{xlim}{Limits of x-axis for the plot.}
	\item{ylim}{Limits of y-axis for the plot.}
	\item{conf.int}{Logical. Whether to display pointwise confidence bands. Default is FALSE.}
	\item{ci.col}{Colour of the confidence bands. Default is \code{col}.}
	\item{ci.lty}{Line type of the confidence bands. Default is 3.}
	\item{legend}{A logical specifying if a legend should be added.}
	\item{legend.pos}{A vector giving the legend's position.
		See \code{\link{legend}} for further details.
	}
	\item{curvlab}{A character or expression vector to appear in the legend.
		Default is the name of the transitions.
	}
	\item{legend.bty}{Box type for the legend. By default no box is drawn.}
	\item{\dots}{Further arguments for plot.}
}

\value{
No value is returned.
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Meira-Machado L., de Uña-Álvarez J., Datta S. (2011). Conditional Transition Probabilities in a non-Markov Illness-death Model. Discussion Papers in Statistics and Operation Research n 11/03. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402-2007). \url{https://depc05.webs.uvigo.es/reports/12_05.pdf}
}

\seealso{
\code{\link{legend}},
\code{\link{plot.default}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2);

# Create survTP object
data(heartTP);
heartTP_obj <- with( heartTP, survTP(time1, event1, Stime, event, age=age) );

# Compute IPCW conditional transition probabilities with confidence band
TPCmsm_obj <- transIPCW(heartTP_obj, s=57, t=310, x=c(0, 15), conf=TRUE, n.boot=100,
method.boot="percentile", method.est=2);

# Build time plots
tr.choice <- dimnames(TPCmsm_obj$est)[[3]];
par.orig <- par( c("mfrow", "cex") );
par( mfrow=c(2,3) );
for ( i in seq_len( length(tr.choice) ) ) {
	plot( TPCmsm_obj, plot.type="t", tr.choice=tr.choice[i], conf.int=TRUE, legend=TRUE,
	main=tr.choice[i], col=seq_len( length(TPCmsm_obj$x) ), lty=1, xlab="", ylab="",
	curvlab=c("Age = 0", "Age = 15") );
}
par(mfrow=c(1, 1), cex=1.2);
title(xlab="Time", ylab="Transition probability", line=3);
par(par.orig);

# Build covariate plots without colors and without confidence band
plot(TPCmsm_obj, plot.type="c", xlab="Age");

# Build covariate plots with colors and without confidence band
plot(TPCmsm_obj, plot.type="c", col=seq_len(5), lty=1, xlab="Age");

# Build covariate plots with confidence band
tr.choice <- dimnames(TPCmsm_obj$est)[[3]];
par.orig <- par( c("mfrow", "cex") );
par( mfrow=c(2,3) );
for ( i in seq_len( length(tr.choice) ) ) {
	plot(TPCmsm_obj, plot.type="c", tr.choice=tr.choice[i], conf.int=TRUE, legend=FALSE,
	main=tr.choice[i], xlab="", ylab="");
}
par(mfrow=c(1, 1), cex=1.2);
title( xlab="Age", ylab=paste("P(", TPCmsm_obj$s, ", ", TPCmsm_obj$t, " | Age)", sep=""), line=3);
par(par.orig);

# Restore the number of threads
setThreadsTP(nth);
}

\keyword{hplot}
\keyword{methods}
\keyword{survival}
