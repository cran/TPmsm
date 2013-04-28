\encoding{UTF-8}
\name{plot.TPmsm}
\alias{plot.TPmsm}
\title{plot method for a TPmsm object}
\description{
	plot method for an object of class 'TPmsm'.
	It draws the estimated transition probabilities in a basic scatterplot.
}
\usage{\S3method{plot}{TPmsm}(x, tr.choice, xlab = "Time", ylab="Transition probability",
col, lty, xlim, ylim, conf.int=FALSE, ci.col, ci.lty,
legend=TRUE, legend.pos, curvlab, legend.bty="n", ...)}

\arguments{
	\item{x}{An object of class 'TPmsm'.}
	\item{tr.choice}{Character vector of the form 'c("from to", "from to")' specifying which transitions should be plotted.
		Default, all the transition probabilities are plotted.}
	\item{xlab}{x-axis label. Default is "Time".}
	\item{ylab}{y-axis label. Default is "Transition probability".}
	\item{col}{Vector of colour. Default is black.}
	\item{lty}{Vector of line type. Default is 1:number of transitions.}
	\item{xlim}{Limits of x-axis for the plot.}
	\item{ylim}{Limits of y-axis for the plot.}
	\item{conf.int}{Logical. Whether to display pointwise confidence bands. Default is FALSE.}
	\item{ci.col}{Colour of the confidence bands. Default is \code{col}.}
	\item{ci.lty}{Line type of the confidence bands. Default is 3.}
	\item{legend}{A logical specifying if a legend should be added.}
	\item{legend.pos}{A vector giving the legend's position. See \code{\link{legend}} for further details.}
	\item{curvlab}{A character or expression vector to appear in the legend. Default is the name of the transitions.}
	\item{legend.bty}{Box type for the legend. By default no box is drawn.}
	\item{\dots}{Further arguments for plot.}
}

\value{
	No value is returned.
}

\author{Artur Agostinho Araújo, Javier Roca-Pardiñas and Luís Meira Machado}

\seealso{
	\code{\link{legend}},
	\code{\link{plot.default}}.
}

\examples{
# set the number of threads
nth <- threadsTP(2)

# Create survTP object
data(bladderTP)
bladderTP_obj <- with(bladderTP, survTP(time1, event1, Stime, event))

# Compute KMW transition probabilities with confidence band
TPmsm_obj <- transKMW(object=bladderTP_obj, s=5, t=59, conf=TRUE, conf.level=0.95,
method.boot="basic", method.est=2)

# Plot all the transitions without confidence band
plot(TPmsm_obj, conf.int=FALSE, col=seq_len(5), lty=1)

# Plot all the transitions with confidence band
tr.choice <- colnames(TPmsm_obj$est)
par.orig <- par( c("mfrow", "cex") )
par( mfrow=c(2,3) )
for ( i in seq_len( length(tr.choice) ) ) {
	plot(TPmsm_obj, tr.choice=tr.choice[i], conf.int=TRUE, legend=FALSE, main=tr.choice[i],
	xlab="", ylab="")
}
par(mfrow=c(1, 1), cex=1.2)
title(xlab="Time", ylab="Transition probability", line=3)
par(par.orig)

# restore the number of threads
threadsTP(nth)
}

\keyword{hplot}
\keyword{methods}
\keyword{survival}
