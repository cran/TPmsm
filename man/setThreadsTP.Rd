\encoding{UTF-8}
\name{setThreadsTP}
\alias{setThreadsTP}
\title{Set the number of threads}

\description{
Specifies the number of threads used by default in parallel sections.
}

\usage{setThreadsTP(num_threads=NULL)}

\arguments{
	\item{num_threads}{the number of threads to use.}
}

\details{
If \code{num_threads} is greater than the number of processors/cores
then the number of processors/cores is used. If package \pkg{TPmsm}
was compiled without OpenMP support then this function returns \code{1}
regardless of the number of processors/cores available.
If \code{num_threads=NULL} the number of threads is not defined.
This is useful when the current number of threads is desired
without defining a new thread number.
}

\value{Invisibly returns the previous number of threads.}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

OpenMP Architecture Review Board, OpenMP Application Program Interface Version 3.0, May 2008, p110. \url{https://www.openmp.org/wp-content/uploads/spec30.pdf}

\dQuote{Runtime Library Routines}, Summary of OpenMP 3.0 C/C++ Syntax, p5. \url{https://www.openmp.org/wp-content/uploads/OpenMP3.0-SummarySpec.pdf}
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\note{
The given thread number is stored in a global variable. This global variable
is then passed to the num_threads clause defined on all parallel sections of
underlying C code. By specifying the number of threads in this way instead of
specifying with a call to omp_set_num_threads we are certain that there is
no interference with the R process. Every time this function is called the
RNG streams are recreated. For more details see \code{\link{setPackageSeedTP}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{environment}
\keyword{utilities}
