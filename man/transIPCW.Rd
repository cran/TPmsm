\name{transIPCW}
\alias{transIPCW}
\title{Inverse probability of censoring weighted transition probabilities}
\description{Provides estimates for the transition probabilities based on inverse probability of censoring weighted estimators, IPCW.}
\usage{transIPCW(object, s, t, x, bw="dpik", window="normal", method.weights="NW",
state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000, conf.level=0.95,
method.boot="percentile", method.est=1, ...)}

\arguments{
	\item{object}{An object of class \code{survTP}.}
	\item{s}{The first time for obtaining estimates for the transition probabilities. If missing, 0 will be used.}
	\item{t}{The second time for obtaining estimates for the transition probabilities. If missing, the maximum of \code{Stime} will be used.}
	\item{x}{Covariate values for obtaining estimates for the conditional transition probabilities. If missing, unconditioned transition probabilities will be computed.}
	\item{bw}{A character string indicating a function to compute a kernel density or nearest-neighbor bandwidth. Defaults to "dpik" from package \pkg{KernSmooth}.}
	\item{window}{A character string specifying the desired kernel or nearest-neighbor window. Defaults to "normal" where the gaussian density kernel will be used.}
	\item{method.weights}{A character string specifying the desired weights method. Possible options are "NW" for the Nadaraya-Watson weights and "LL" for local linear weights. Defaults to "NW".}
	\item{state.names}{A vector of characters giving the state names.}
	\item{conf}{Provides pointwise confidence bands. Defaults to \code{FALSE}.}
	\item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
	\item{conf.level}{Level of confidence. Defaults to 0.95 (correponding to 95\%).}
	\item{method.boot}{The method used to compute bootstrap confidence bands. Possible options are "percentile" and "basic". Defaults to "percentile".}
	\item{method.est}{The method used to compute the estimate. Possible options are 1 or 2.}
	\item{\dots}{Further arguments. Typically these arguments are passed to the function specified by argument \code{bw}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
	If argument \code{x} is missing or if argument \code{object} doesn't contain a covariate,
	an object of class \code{TPmsm} is returned. There are methods for \code{contour}, \code{image}, \code{print} and \code{plot}.
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

	If argument \code{x} is specified and argument \code{object} contains a covariate,
	an object of class \code{TPCmsm} is returned. There are methods for \code{print} and \code{plot}.
	\code{TPCmsm} objects are implemented as a list with elements:
	\item{method}{A string indicating the type of estimator used in the computation.}
	\item{est}{A 3 dimentional array with transition probability estimates. The first dimension being the event times, the second the covariate values and the last one the 5 possible transitions.}
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

\author{Artur Agostinho Araujo, Javier Roca-Pardinas and Luis Meira Machado}

\references{
	Luis Meira-Machado, Jacobo de Una-Alvarez and Somnath Datta. Conditional Transition Probabilities in a non-Markov Illness-death Model. Discussion Papers in Statistics and Operation Research n 11/03, 2011. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402 - 2007). This file can be downloaded from: \url{http://webs.uvigo.es/depc05/reports/12_05.pdf}

	Meira Machado LF, de Una-Alvarez J, Cadarso-Suarez C (2006). "Nonparametric estimation of transition probabilities in a non-Markov illness-death model", \emph{Lifetime Data Anal}, \bold{12}, 325-344.

	Davison, A.C. and Hinkley, D.V. (1997) "Bootstrap Methods and Their Application", Chapter 5. Cambridge University Press.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{transAJ}},
	\code{\link{transKMPW}},
	\code{\link{transKMW}},
	\code{\link{transLIN}},
	\code{\link{transLS}},
	\code{\link{transPAJ}}.
}

\examples{
# Create survTP object with age as covariate
data(heartTP)
heartTP_obj <- with(heartTP, survTP(time1, event1, Stime, event, age))

# Compute unconditioned transition probabilities
transIPCW(object=heartTP_obj, s=33, t=412)

# Compute unconditioned transition probabilities with confidence band
transIPCW(object=heartTP_obj, s=33, t=412, conf=TRUE, conf.level=0.9,
method.boot="basic", method.est=2)

# Compute conditional transition probabilities
transIPCW(object=heartTP_obj, s=33, t=412, x=0)

# Compute conditional transition probabilities with confidence band
transIPCW(object=heartTP_obj, s=33, t=412, x=0, conf=TRUE, conf.level=0.95,
method.boot="percentile", method.est=2)
}
