\encoding{UTF-8}
\name{bladderTP}
\alias{bladderTP}
\docType{data}
\title{Bladder cancer recurrences}

\description{
This contains the bladder cancer recurrences data in a different format.
In this study, patients had superficial bladder tumors that were removed by transurethral resection.
Many patients had multiple recurrences (up to a maximum of 9) of tumors during the study,
and new tumors were removed at each visit. Only the first two recurrence times (in months) are considered.
}

\usage{data(bladderTP)}

\format{
	A data frame with 85 observations on the following 4 variables.
	\describe{
		\item{\code{time1}}{Time until first recurrence/censoring time.}
		\item{\code{event1}}{First recurrence indicator.}
		\item{\code{Stime}}{Time until second recurrence/censoring time, i.e. the total time of the process.}
		\item{\code{event}}{First or second recurrence indicator, i.e. the censoring indicator of the total time.}
	}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
Wei L. J., Lin D. Y., Weissfeld L. (1989) Regression analysis of multivariate incomplete failure time data by modeling marginal distributions. \emph{Journal of the American Statistical Association} \bold{84}(408), 1065--1073. \doi{10.2307/2290084}
}

%\examples{
%data(bladderTP)
%## maybe str(bladderTP) ; plot(bladderTP) ...
%}

\keyword{datasets}
