\encoding{UTF-8}
\name{threadsTP}
\alias{threadsTP}
\title{Set the number of threads}

\description{
  Specifies the number of threads used by default in parallel sections.
}

\usage{threadsTP(num_threads=NULL)}

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
  OpenMP Architecture Review Board, OpenMP Application Program Interface Version 3.0, May 2008, p110 (\url{http://www.openmp.org/mp-documents/spec30.pdf})\cr
	\dQuote{Runtime Library Routines}, Summary of OpenMP 3.0 C/C++ Syntax, p5 (\url{http://www.openmp.org/mp-documents/OpenMP3.0-SummarySpec.pdf})\cr
}

\author{Artur Agostinho Araújo, Javier Roca-Pardiñas and Luís Meira Machado}

\note{
  The given thread number is stored in a global variable. This global variable
  is then passed to the num_threads clause defined on all parallel sections of
  underlying C code. By specifying the number of threads in this way instead of
  specifying with a call to omp_set_num_threads we are certain that there is
  no interference with the R process.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}

\examples{
# set the number of threads
nth <- threadsTP(2)

# restore the number of threads
threadsTP(nth)
}

\keyword{environment}
\keyword{utilities}
