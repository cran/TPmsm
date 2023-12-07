\encoding{UTF-8}
\name{TPmsmOut}
\alias{TPmsmOut}
\title{Convert a data.frame in the \pkg{TPmsm} format to other formats}

\description{Converts a data.frame in the \pkg{TPmsm} format to formats supported by external packages.}

\usage{TPmsmOut(data, names, package="p3state.msm")}

\arguments{
	\item{data}{A data.frame in the \pkg{TPmsm} format.}
	\item{names}{A character vector of lenght 4,
		indicating the variable names equivalent to variable names
		\dQuote{time1}, \dQuote{event1}, \dQuote{Stime}, \dQuote{event}
		in the \pkg{TPmsm} format, in this order.
	}
	\item{package}{The format of the data.frame.
		Possible options are \dQuote{p3state.msm} and \dQuote{etm}.
		Defaults to \dQuote{p3state.msm}.
	}
}

\value{A data.frame in the format specified by argument \code{package}.}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\references{
Allignol A., Schumacher M., Beyersmann J. (2011). Empirical Transition Matrix of Multi-State Models: The etm Package. \emph{Journal of Statistical Software}, \bold{38}(4), 1-15. \doi{10.18637/jss.v038.i04}

Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Meira-Machado L., Roca-Pardiñas J. (2011). p3state.msm: Analyzing Survival Data from an Illness-Death Model. \emph{Journal of Statistical Software}, \bold{38}(3), 1-18. \doi{10.18637/jss.v038.i03}
}

\seealso{
\code{\link{as.data.frame.survTP}},
\code{\link{survTP}}.
}

\examples{
data(heartTP);
heartP3 <- TPmsmOut( heartTP, c("time1", "event1", "Stime", "event") );
head(heartP3);
}

\keyword{manip}
\keyword{survival}
