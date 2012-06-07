\name{transKMPW}
\alias{transKMPW}
\title{Presmoothed Kaplan-Meier weighted transition probabilities}
\description{Provides estimates for the transition probabilities based on presmoothed Kaplan-Meier weighted estimators, KMPW.}
\usage{transKMPW(object, s, t, state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000,
conf.level=0.95, method.boot="percentile", method.est=3)}

\arguments{
	\item{object}{An object of class \code{survTP}.}
	\item{s}{The first time for obtaining estimates for the transition probabilities. If missing, 0 will be used.}
	\item{t}{The second time for obtaining estimates for the transition probabilities. If missing, the maximum of \code{Stime} will be used.}
	\item{state.names}{A vector of characters giving the state names.}
	\item{conf}{Provides pointwise confidence bands. Defaults to \code{FALSE}.}
	\item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
	\item{conf.level}{Level of confidence. Defaults to 0.95 (correponding to 95\%).}
	\item{method.boot}{The method used to compute bootstrap confidence bands. Possible options are "percentile" and "basic". Defaults to "percentile".}
	\item{method.est}{The method used to compute the estimate. Possible options are 1, 2, 3 or 4.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
	An object of class \code{TPmsm}. There are methods for \code{contour}, \code{image}, \code{print} and \code{plot}.
	\code{TPmsm} objects are implemented as a list with elements:
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

\author{Artur Agostinho Araujo, Javier Roca-Pardinas and Luis Meira Machado}

\references{
	Amorim AP, de Una-Alvarez J, Meira Machado LF (2011). "Presmoothing the transition probabilities in the illness-death model", \emph{Statistics and Probability Letters}, \bold{81}, 797-806.

	Davison, A.C. and Hinkley, D.V. (1997) "Bootstrap Methods and Their Application", Chapter 5. Cambridge University Press.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{transAJ}},
	\code{\link{transIPCW}},
	\code{\link{transKMW}},
	\code{\link{transLIN}},
	\code{\link{transLS}},
	\code{\link{transPAJ}}.
}

\examples{
# Create survTP object
data(heartTP)
heartTP_obj <- with(heartTP, survTP(time1, event1, Stime, event))

# Compute transition probabilities
transKMPW(object=heartTP_obj, s=33, t=412)

# Compute transition probabilities with confidence band
transKMPW(object=heartTP_obj, s=33, t=412, conf=TRUE, conf.level=0.9,
method.boot="percentile", method.est=4)
}
