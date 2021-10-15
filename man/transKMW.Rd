\encoding{UTF-8}
\name{transKMW}
\alias{transKMW}
\title{Kaplan-Meier weighted transition probabilities}

\description{Provides estimates for the transition probabilities based on Kaplan-Meier weighted estimators, KMW.}

\usage{transKMW(object, s, t, state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000,
conf.level=0.95, method.boot="percentile", method.est=3)}

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
	\item{method.est}{The method used to compute the estimate. Possible options are 1, 2, 3 or 4.}
}

\details{
If \code{method.est=1} then \eqn{p_{11}(s,t)}, \eqn{p_{12}(s,t)} and \eqn{p_{22}(s,t)} are estimated according to the following expressions:\cr\cr
\eqn{p_{11}(s,t)=\frac{1-P(Z \leq t)}{1-P(Z \leq s)}},\cr\cr
\eqn{p_{12}(s,t)=\frac{P(Z \leq t)-P(Z \leq s)-P(s<Z \leq t, T \leq t)}{1-P(Z \leq s)}},\cr\cr
\eqn{p_{22}(s,t) =\frac{P(Z \leq s)-P(Z \leq s,T \leq t)}{P(Z \leq s)-P(T \leq s)}}.\cr\cr
Then, \eqn{p_{13}(s,t)=1-p_{11}(s,t)-p_{12}(s,t)} and \eqn{p_{23}(s,t)=1-p_{22}(s,t)}.\cr\cr
If \code{method.est=2} then \eqn{p_{11}(s,t)}, \eqn{p_{12}(s,t)} and \eqn{p_{22}(s,t)} are estimated according to the following expressions:\cr\cr
\eqn{p_{11}(s,t)=\frac{P(Z>t)}{P(Z>s)}},\cr\cr
\eqn{p_{12}(s,t)=\frac{P(s<Z \leq t,T>t)}{P(Z>s)}},\cr\cr
\eqn{p_{22}(s,t) =\frac{P(Z \leq s,T>t)}{P(Z \leq s,T>s)}}.\cr\cr
Then, \eqn{p_{13}(s,t)=1-p_{11}(s,t)-p_{12}(s,t)} and \eqn{p_{23}(s,t)=1-p_{22}(s,t)}.\cr\cr
If \code{method.est=3} then \eqn{p_{11}(s,t)}, \eqn{p_{13}(s,t)} and \eqn{p_{23}(s,t)} are estimated according to the following expressions:\cr\cr
\eqn{p_{11}(s,t)=\frac{1-P(Z \leq t)}{1-P(Z \leq s)}},\cr\cr
\eqn{p_{13}(s,t)=\frac{P(Z>s,T \leq t)}{1-P(Z \leq s)}},\cr\cr
\eqn{p_{23}(s,t) =\frac{P(Z \leq s,s<T \leq t)}{P(Z \leq s)-P(T \leq s)}}.\cr\cr
Then, \eqn{p_{12}(s,t)=1-p_{11}(s,t)-p_{13}(s,t)} and \eqn{p_{22}(s,t)=1-p_{23}(s,t)}.\cr\cr
If \code{method.est=4} then \eqn{p_{11}(s,t)}, \eqn{p_{13}(s,t)} and \eqn{p_{23}(s,t)} are estimated according to the following expressions:\cr\cr
\eqn{p_{11}(s,t)=\frac{P(Z>t)}{P(Z>s)}},\cr\cr
\eqn{p_{13}(s,t)=\frac{P(Z>s,T \leq t)}{P(Z>s)}},\cr\cr
\eqn{p_{23}(s,t) =\frac{P(Z \leq s,s<T \leq t)}{P(Z \leq s,T>s)}}.\cr\cr
Then, \eqn{p_{12}(s,t)=1-p_{11}(s,t)-p_{13}(s,t)} and \eqn{p_{22}(s,t)=1-p_{23}(s,t)}.
}

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
\code{\link{transLIN}},
\code{\link{transLS}},
\code{\link{transPAJ}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Create survTP object
data(heartTP)
heartTP_obj <- with(heartTP, survTP(time1, event1, Stime, event))

# Compute transition probabilities
transKMW(object=heartTP_obj, s=33, t=412)

# Compute transition probabilities with confidence band
transKMW(object=heartTP_obj, s=33, t=412, conf=TRUE, conf.level=0.9,
method.boot="basic", method.est=2)

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{dplot}
\keyword{nonparametric}
\keyword{survival}
