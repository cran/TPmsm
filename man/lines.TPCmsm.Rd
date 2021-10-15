\encoding{UTF-8}
\name{lines.TPCmsm}
\alias{lines.TPCmsm}
\title{lines method for a TPCmsm object}

\description{lines method for an object of class \sQuote{TPCmsm}.}

\usage{\S3method{lines}{TPCmsm}(x, plot.type="t", tr.choice, col, lty, conf.int=FALSE,
ci.col, ci.lty, legend=FALSE, legend.pos, curvlab, legend.bty="n", ...)}

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
	\item{col}{Vector of colour.}
	\item{lty}{Vector of line type. Default is 1:number of transitions.}
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
	\item{\dots}{Further arguments for lines.}
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
\code{\link{lines}},
\code{\link{plot.default}},
\code{\link{plot.TPCmsm}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Create survTP object
data(heartTP)
heartTP_obj <- with( heartTP, survTP(time1, event1, Stime, event, age=age) )

# Compute IPCW1 conditional transition probabilities without confidence band
TPC_IPCW1 <- transIPCW(heartTP_obj, s=57, t=310, x=15, conf=FALSE, method.est=1)

# Compute IPCW2 conditional transition probabilities without confidence band
TPC_IPCW2 <- transIPCW(heartTP_obj, s=57, t=310, x=15, conf=FALSE, method.est=2)

# Compute LIN conditional transition probabilities without confidence band
TPC_LIN <- transLIN(heartTP_obj, s=57, t=310, x=15, conf=FALSE)

# Build covariate plots
tr.choice <- dimnames(TPC_LIN$est)[[3]]
par.orig <- par( c("mfrow", "cex") )
par( mfrow=c(2,3) )
for ( i in seq_len( length(tr.choice) ) ) {
	plot(TPC_IPCW1, plot.type="c", tr.choice=tr.choice[i], legend=FALSE,
	main=tr.choice[i], col=1, lty=1, xlab="", ylab="")
	lines(TPC_IPCW2, plot.type="c", tr.choice=tr.choice[i], legend=FALSE, col=2, lty=1)
	lines(TPC_LIN, plot.type="c", tr.choice=tr.choice[i], legend=FALSE, col=3, lty=1)
}
plot.new()
legend(x="center", legend=c("IPCW1", "IPCW2", "LIN"), col=1:3, lty=1, bty="n", cex=1.5)
par(mfrow=c(1, 1), cex=1.2)
title(xlab="Age", ylab="Transition probability", line=3)
par(par.orig)

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{aplot}
\keyword{methods}
\keyword{survival}
