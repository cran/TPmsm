\name{survTP}
\alias{survTP}
\alias{is.survTP}
\title{Create a survTP object}
\description{Creates a survTP object, usually used as input to other functions.}
\usage{
survTP(time1, event1, Stime, event, covariate=NULL)
is.survTP(object)
}

\arguments{
\item{time1}{Time of the transition into state 2, state 3 or censoring time.}
\item{event1}{Indicator of transition into state 2 or state 3; 0 if the transition time is censored and 1 otherwise.}
\item{Stime}{The total time of the process.}
\item{event}{Censoring indicator of the total time of the process; 0 if the total time is censored and 1 otherwise.}
\item{covariate}{Any covariate variable.}
\item{object}{Any R object.}
}

%\details{
%  Typical usages are \preformatted{survTP(time1, event1, Stime, event)}
%}

\value{
An object of class \code{survTP}.

\code{survTP} objects are implemented as a single element list \item{data}{a data.frame
with \code{time1}, \code{event1}, \code{Stime}, \code{event} and \code{covariate} as columns.}

In the case of \code{is.survTP}, a logical value \code{TRUE} if \code{object}
inherits from class \code{survTP}, otherwise \code{FALSE}.
}

\author{Artur Agostinho Araujo, Javier Roca-Pardinas and Luis Meira Machado}

\examples{
data(bladderTP)
bladderTP_obj <- with(bladderTP, survTP(time1, event1, Stime, event))
#or
bladderTP_obj <- survTP(bladderTP$time1, bladderTP$event1, bladderTP$Stime,
bladderTP$event)
data(heartTP)
heartTP_obj <- with(heartTP, survTP(time1, event1, Stime, event, age))
#or
heartTP_obj <- survTP(heartTP$time1, heartTP$event1, heartTP$Stime,
heartTP$event, heartTP$age)
}
