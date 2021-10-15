\encoding{UTF-8}
\name{transPAJ}
\alias{transPAJ}
\title{Presmoothed Aalen-Johansen transition probabilities}

\description{Provides estimates for the transition probabilities based on the presmoothed Aalen-Johansen estimator, PAJ.}

\usage{transPAJ(object, s, t, state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000,
conf.level=0.95, method.boot="percentile")}

\arguments{
	\item{object}{An object of class \sQuote{survTP}.}
	\item{s}{The first time for obtaining estimates for the transition probabilities.
		If missing, 0 will be used.
	}
	\item{t}{The second time for obtaining estimates for the transition probabilities.
		If missing, the maximum of \code{Stime} will be used.
	}
	\item{state.names}{A vector of characters giving the state names.}
	\item{conf}{Provides pointwise confidence bands. Defaults to \code{FALSE}.}
	\item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
	\item{conf.level}{Level of confidence. Defaults to 0.95 (corresponding to 95\%).}
	\item{method.boot}{The method used to compute bootstrap confidence bands.
		Possible options are \dQuote{percentile} and \dQuote{basic}.
		Defaults to \dQuote{percentile}.
	}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
An object of class \sQuote{TPmsm}. There are methods for \code{contour}, \code{image}, \code{print} and \code{plot}.
\sQuote{TPmsm} objects are implemented as a list with elements:
\item{method}{A string indicating the type of estimator used in the computation.}
\item{est}{A matrix with transition probability estimates. The rows being the event times and the columns the 5 possible transitions.}
\item{inf}{A matrix with the lower transition probabilities of the confidence band. The rows being the event times and the columns the 5 possible transitions.}
\item{sup}{A matrix with the upper transition probabilities of the confidence band. The rows being the event times and the columns the 5 possible transitions.}
\item{time}{Vector of times where the transition probabilities are computed.}
\item{s}{Start of the time interval.}
\item{t}{End of the time interval.}
\item{h}{The bandwidth used. If the estimator doesn't require a bandwidth, it's set to \code{NULL}.}
\item{state.names}{A vector of characters giving the states names.}
\item{n.boot}{Number of bootstrap samples used in the computation of the confidence band.}
\item{conf.level}{Level of confidence used to compute the confidence band.}
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Moreira A., de Uña-Álvarez J. and Meira-Machado L. (2011). Presmoothing the Aalen-Johansen estimator of transition probabilities. Discussion Papers in Statistics and Operation Research n 11/03. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402-2007). \url{https://depc05.webs.uvigo.es/reports/11_03.pdf}

Davison A. C., Hinkley D. V. (1997). \emph{Bootstrap Methods and their Application}, Chapter 5, Cambridge University Press.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{transAJ}},
\code{\link{transIPCW}},
\code{\link{transKMPW}},
\code{\link{transKMW}},
\code{\link{transLIN}},
\code{\link{transLS}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Create survTP object
data(heartTP)
heartTP_obj <- with(heartTP, survTP(time1, event1, Stime, event))

# Compute transition probabilities
transPAJ(object=heartTP_obj, s=33, t=412)

# Compute transition probabilities with confidence band
transPAJ(object=heartTP_obj, s=33, t=412, conf=TRUE, conf.level=0.9,
method.boot="percentile")

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{dplot}
\keyword{models}
\keyword{regression}
\keyword{survival}
