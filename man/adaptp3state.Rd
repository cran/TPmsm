\name{adaptp3state}
\alias{adaptp3state}
\title{Convert a data.frame in the \pkg{TPmsm} format to the \pkg{p3state.msm} format}
\description{Converts a data.frame in the \pkg{TPmsm} format to the \pkg{p3state.msm} format.}
\usage{adaptp3state(data, varnames)}

\arguments{
	\item{data}{A data.frame in the \pkg{TPmsm} format.}
	\item{varnames}{
		A character vector of lenght 4,
		indicating the variable names equivalent to variable names
		"time1", "event1", "Stime", "event" in the \pkg{TPmsm} format,
		in this order.
	}
}

%\details{
%  Typical usages are \preformatted{adaptp3state(data, varnames)}
%}

\value{A data.frame in the \pkg{p3state.msm} format.}

\author{Artur Agostinho Araujo, Javier Roca-Pardinas and Luis Meira Machado}

\references{
	Meira-Machado L, Roca-Pardinas J (2011). "p3state.msm: Analyzing Survival Data from an Illness-Death Model." Journal of Statistical Software, 38(3): 1-18.
}

\examples{
data(heartTP)
heartP3 <- adaptp3state(heartTP, c("time1", "event1", "Stime", "event"))
}
