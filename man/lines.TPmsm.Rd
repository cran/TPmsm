\encoding{UTF-8}
\name{lines.TPmsm}
\alias{lines.TPmsm}
\title{lines method for a TPmsm object}

\description{lines method for an object of class \sQuote{TPmsm}.}

\usage{\S3method{lines}{TPmsm}(x, tr.choice, col, lty, conf.int=FALSE, ci.col, ci.lty,
legend=FALSE, legend.pos, curvlab, legend.bty="n", ...)}

\arguments{
	\item{x}{An object of class \sQuote{TPmsm}.}
	\item{tr.choice}{Character vector of the form \sQuote{c(\dQuote{from to}, \dQuote{from to})}
		specifying which transitions should be plotted.
		Default, all the transition probabilities are plotted.
	}
	\item{col}{Vector of colour. Default is black.}
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
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014) TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}
}

\seealso{
\code{\link{legend}},
\code{\link{lines}},
\code{\link{plot.default}},
\code{\link{plot.TPmsm}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Create survTP object
data(bladderTP)
bladderTP_obj <- with( bladderTP, survTP(time1, event1, Stime, event) )

# Compute transition probabilities without confidence band
KMW <- transKMW(object=bladderTP_obj, s=5, t=59, conf=FALSE, method.est=1)
KMPW <- transKMPW(object=bladderTP_obj, s=5, t=59, conf=FALSE, method.est=1)
AJ <- transAJ(object=bladderTP_obj, s=5, t=59, conf=FALSE)
PAJ <- transPAJ(object=bladderTP_obj, s=5, t=59, conf=FALSE)
LIN <- transLIN(object=bladderTP_obj, s=5, t=59, conf=FALSE)
LS <- transLS(object=bladderTP_obj, s=5, t=59, h=c(0.25, 2.5),
nh=25, ncv=50, conf=FALSE)

# Plot '1 2' KMW transition probability estimate
par( mfrow=c(1, 1) )
plot(KMW, tr.choice="1 2", ylab="P12(5, Time)", xlab="Time",
col=1, lty=1, legend=FALSE)

# Add other '1 2' transition probability estimates
lines(KMPW, tr.choice="1 2", col=2, lty=1)
lines(AJ, tr.choice="1 2", col=3, lty=1)
lines(PAJ, tr.choice="1 2", col=4, lty=1)
lines(LIN, tr.choice="1 2", col=5, lty=1)
lines(LS, tr.choice="1 2", col=6, lty=1)

# Add legend
legend(x="topleft", legend=c("KMW", "KMPW", "AJ", "PAJ", "LIN", "LS"),
col=1:6, lty=1, bty="n")

# Plot all the transitions
tr.choice <- colnames(KMW$est)
par.orig <- par( c("mfrow", "cex") )
par( mfrow=c(2, 3) )
for ( i in seq_len( length(tr.choice) ) ) {
	plot(KMW, tr.choice=tr.choice[i], col=1, lty=1, legend=FALSE,
	main=tr.choice[i], xlab="", ylab="")
	lines(KMPW, tr.choice=tr.choice[i], col=2, lty=1)
	lines(AJ, tr.choice=tr.choice[i], col=3, lty=1)
	lines(PAJ, tr.choice=tr.choice[i], col=4, lty=1)
	lines(LIN, tr.choice=tr.choice[i], col=5, lty=1)
	lines(LS, tr.choice=tr.choice[i], col=6, lty=1)
}
plot.new()
legend(x="center", legend=c("KMW", "KMPW", "AJ", "PAJ", "LIN", "LS"),
col=1:6, lty=1, bty="n", cex=1.5)
par(mfrow=c(1, 1), cex=1.2)
title(xlab="Time", ylab="Transition probability", line=3)
par(par.orig)

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{aplot}
\keyword{methods}
\keyword{survival}
