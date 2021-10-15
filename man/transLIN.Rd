\encoding{UTF-8}
\name{transLIN}
\alias{transLIN}
\title{LIN based transition probabilities}

\description{Provides estimates for the transition probabilities based on LIN estimators, LIN.}

\usage{transLIN(object, s, t, x, bw="dpik", window="normal", method.weights="NW",
state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000, conf.level=0.95,
method.boot="percentile", ...)}

\arguments{
	\item{object}{An object of class \sQuote{survTP}.}
	\item{s}{The first time for obtaining estimates for the transition probabilities.
		If missing, 0 will be used.
	}
	\item{t}{The second time for obtaining estimates for the transition probabilities.
		If missing, the maximum of \code{Stime} will be used.
	}
	\item{x}{Covariate values for obtaining estimates for the conditional transition probabilities.
		If missing, unconditioned transition probabilities will be computed.
	}
	\item{bw}{A character string indicating a function to compute a kernel density bandwidth.
		Defaults to \dQuote{dpik} from package \pkg{KernSmooth}.
		Alternatively a single numeric value can be specified.
	}
	\item{window}{A character string specifying the desired kernel.
		See details below for possible options.
		Defaults to \dQuote{normal} where the gaussian density kernel will be used.
	}
	\item{method.weights}{A character string specifying the desired weights method.
		Possible options are \dQuote{NW} for the Nadaraya-Watson weights
		and \dQuote{LL} for local linear weights. Defaults to \dQuote{NW}.
	}
	\item{state.names}{A vector of characters giving the state names.}
	\item{conf}{Provides pointwise confidence bands. Defaults to \code{FALSE}.}
	\item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
	\item{conf.level}{Level of confidence. Defaults to 0.95 (corresponding to 95\%).}
	\item{method.boot}{The method used to compute bootstrap confidence bands.
		Possible options are \dQuote{percentile} and \dQuote{basic}.
		Defaults to \dQuote{percentile}.
	}
	\item{\dots}{Further arguments.
		Typically these arguments are passed to the function specified by argument \code{bw}.
	}
}

\details{
If \code{bw="dpik"} then possible options for argument \code{window} are \dQuote{normal}, \dQuote{box}, \dQuote{epanech}, \dQuote{biweight} or \dQuote{triweight}.
When argument \code{bw} is numeric then argument \code{window} accepts the same options as when \code{bw="dpik"} plus one of \dQuote{tricube}, \dQuote{triangular} or \dQuote{cosine}.
}

\value{
If argument \code{x} is missing or if argument \code{object} doesn't contain a covariate,
an object of class \sQuote{TPmsm} is returned. There are methods for \code{contour}, \code{image}, \code{print} and \code{plot}.
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

If argument \code{x} is specified and argument \code{object} contains a covariate,
an object of class \sQuote{TPCmsm} is returned. There are methods for \code{print} and \code{plot}.
\sQuote{TPCmsm} objects are implemented as a list with elements:
\item{method}{A string indicating the type of estimator used in the computation.}
\item{est}{A 3 dimensional array with transition probability estimates. The first dimension being the event times, the second the covariate values and the last one the 5 possible transitions.}
\item{inf}{A 3 dimensional array with the lower transition probabilities of the confidence band. The first dimension being the event times, the second the covariate values and the last one the 5 possible transitions.}
\item{sup}{A 3 dimensional array with the upper transition probabilities of the confidence band. The first dimension being the event times, the second the covariate values and the last one the 5 possible transitions.}
\item{time}{Vector of times where the transition probabilities are computed.}
\item{covariate}{Vector of covariate values where the conditional transition probabilities are computed.}
\item{s}{Start of the time interval.}
\item{t}{End of the time interval.}
\item{x}{Additional covariate values where the conditional transition probabilities are computed, which may or may not be present in the sample.}
\item{h}{The bandwidth used.}
\item{state.names}{A vector of characters giving the states names.}
\item{n.boot}{Number of bootstrap samples used in the computation of the confidence band.}
\item{conf.level}{Level of confidence used to compute the confidence band.}
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Meira-Machado L., de Uña-Álvarez J., Datta S. (2011). Conditional Transition Probabilities in a non-Markov Illness-death Model. Discussion Papers in Statistics and Operation Research n 11/03. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402-2007). \url{https://depc05.webs.uvigo.es/reports/12_05.pdf}

Meira Machado L. F., de Uña-Álvarez J., Cadarso-Suárez C. (2006). Nonparametric estimation of transition probabilities in a non-Markov illness-death model. \emph{Lifetime Data Anal}, \bold{12}(3), 325-344. \doi{10.1007/s10985-006-9009-x}

Davison, A. C., Hinkley, D. V. (1997). \emph{Bootstrap Methods and their Application}, Chapter 5, Cambridge University Press.
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
\code{\link{transLS}},
\code{\link{transPAJ}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Create survTP object with age as covariate
data(heartTP)
heartTP_obj <- with(heartTP, survTP(time1, event1, Stime, event, age=age))

# Compute unconditioned transition probabilities
transLIN(object=heartTP_obj, s=33, t=412)

# Compute unconditioned transition probabilities with confidence band
transLIN(object=heartTP_obj, s=33, t=412, conf=TRUE, conf.level=0.9,
method.boot="basic")

# Compute conditional transition probabilities
transLIN(object=heartTP_obj, s=33, t=412, x=0)

# Compute conditional transition probabilities with confidence band
transLIN(object=heartTP_obj, s=33, t=412, x=0, conf=TRUE, conf.level=0.95,
n.boot=100, method.boot="percentile")

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{smooth}
\keyword{survival}
