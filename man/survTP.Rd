\encoding{UTF-8}
\name{survTP}
\alias{survTP}
\alias{is.survTP}
\title{Create a survTP object}

\description{Creates a \sQuote{survTP} object, usually used as input to other functions.}

\usage{
	survTP(time1, event1, Stime, event, ...)
	is.survTP(x)
}

\arguments{
	\item{time1}{Time of the transition into state 2, state 3 or censoring time.}
	\item{event1}{Indicator of transition into state 2 or state 3; 0 if the transition time is censored and 1 otherwise.}
	\item{Stime}{The total time of the process.}
	\item{event}{Censoring indicator of the total time of the process; 0 if the total time is censored and 1 otherwise.}
	\item{\dots}{Any number of covariates can be specified.}
	\item{x}{Any R object.}
}

%\details{
%  Typical usages are \preformatted{survTP(time1, event1, Stime, event, ...)}
%}

\value{
An object of class \sQuote{survTP}.

\sQuote{survTP} objects are implemented as a single element list \item{data}{a data.frame
with \code{time1}, \code{event1}, \code{Stime}, \code{event} and covariates as columns.}

In the case of \code{is.survTP}, a logical value \code{TRUE} if \code{x}
inherits from class \sQuote{survTP}, otherwise \code{FALSE}.
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}
}

\examples{
data(bladderTP);
bladderTP_obj <- with( bladderTP, survTP(time1, event1, Stime, event) );
#or
bladderTP_obj <- survTP(bladderTP$time1, bladderTP$event1, bladderTP$Stime,
bladderTP$event);
data(heartTP);
heartTP_obj <- with( heartTP, survTP(time1, event1, Stime, event, age=age) );
#or
heartTP_obj <- survTP(heartTP$time1, heartTP$event1, heartTP$Stime,
heartTP$event, age=heartTP$age);
}

\keyword{multivariate}
\keyword{survival}
