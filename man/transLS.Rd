\encoding{UTF-8}
\name{transLS}
\alias{transLS}
\title{Location-Scale transition probabilities}
\description{Provides estimates for the transition probabilities based on the Location-Scale estimator, LS.}
\usage{transLS(object, s, t, h, nh=40, ncv=10, window="normal", state.names=c("1", "2", "3"),
conf=FALSE, n.boot=1000, conf.level=0.95, method.boot="percentile", boot.cv=FALSE,
cv.full=TRUE)}

\arguments{
	\item{object}{An object of class \code{survTP}.}
	\item{s}{The first time for obtaining estimates for the transition probabilities. If missing, 0 will be used.}
	\item{t}{The second time for obtaining estimates for the transition probabilities. If missing, the maximum of \code{Stime} will be used.}
	\item{h}{A vector with 1 up to 4 values, indicating the minimum and maximum bandwidths to test by cross-validation.}
	\item{nh}{The number of bandwidth values to test by cross-validation. Defaults to 40.}
	\item{ncv}{The number of cross-validation samples. Defaults to 10.}
	\item{window}{A character string specifying the desired kernel.
		Possible options are "normal", "epanech", "biweight", "triweight", "box", "tricube", "triangular" or "cosine".
		Defaults to "normal" where the gaussian density kernel will be used.}
	\item{state.names}{A vector of characters giving the state names.}
	\item{conf}{Provides pointwise confidence bands. Defaults to \code{FALSE}.}
	\item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
	\item{conf.level}{Level of confidence. Defaults to 0.95 (corresponding to 95\%).}
	\item{method.boot}{The method used to compute bootstrap confidence bands. Possible options are "percentile" and "basic". Defaults to "percentile".}
	\item{boot.cv}{If \code{TRUE} the bandwidth is computed by cross-validation for each bootstrap sample.
		If \code{FALSE} the bandwidth used to compute the estimates is used to compute each bootstrap estimate. Defaults to \code{FALSE}.}
	\item{cv.full}{If \code{TRUE} the bandwidth is computed by cross-validation for both the location and scale functions.
		If \code{FALSE} the bandwidth is computed by cross-validation only for the location function.
		And the bandwidth for the scale function is taken to be equal to the location one. Defaults to \code{TRUE}.}
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

\author{Artur Agostinho Araújo, Javier Roca-Pardiñas and Luís Meira Machado}

\references{
	Meira-Machado L., Roca-Pardiñas J., Van Keilegom I., and Cadarso-Suárez C. (2013). Bandwidth Selection for the Estimation of Transition Probabilities in the Location-Scale Progressive Three-State Model. \emph{Computational Statistics}, \bold{28}, 2185-2210.

	Meira-Machado L., Roca-Pardiñas J., Van Keilegom I. and Cadarso-Suárez C. Estimation of transition probabilities in a non-Markov model with successive survival times. Discussion paper 2010. This file can be downloaded from: \url{http://sites.uclouvain.be/IAP-Stat-Phase-V-VI/ISBApub/dp2010/DP1053.pdf}

	Van Keilegom I., de Uña-Álvarez J. and Meira-Machado L. (2011) Nonparametric location-scale models for successive survival times under dependent censoring. \emph{Journal of Statistical Planning and Inference}, \bold{141}, 1118-1131.

	Davison, A. C. and Hinkley, D. V. (1997) Bootstrap Methods and their Application, Chapter 5. Cambridge University Press.
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
	\code{\link{transPAJ}}.
}

\examples{
# set the number of threads
nth <- threadsTP(2)

# Create survTP object
data(bladderTP)
bladderTP_obj <- with(bladderTP, survTP(time1, event1, Stime, event))

# Compute transition probabilities
LS0 <- transLS(object=bladderTP_obj, s=5, t=59, h=c(0.25, 2.5), nh=25, ncv=50, conf=FALSE)
print(LS0)

# Compute transition probabilities with confidence band
h <- with( LS0, c( rep(h[1], 2), rep(h[2], 2) ) )
transLS(object=bladderTP_obj, s=5, t=59, h=h, conf=TRUE,
conf.level=0.95, method.boot="percentile", boot.cv=FALSE)

# restore the number of threads
threadsTP(nth)
}

\keyword{dplot}
\keyword{nonparametric}
\keyword{smooth}
\keyword{survival}
