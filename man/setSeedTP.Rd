\encoding{UTF-8}
\name{setSeedTP}
\alias{setSeedTP}
\title{Save and restore RNG stream seeds}

\description{
The random number generator (RNG) with multiple independent streams developed by
L'Ecuyer et al. (2002) is used for parallel computation of uniform pseudorandom numbers.
Package \pkg{TPmsm} makes extensive use of uniform pseudorandom numbers,
particularly for the bootstrapping statistical techniques and for
the generation of univariate and multivariate pseudorandom data.
This function permits saving and restoring the seed of each individual RNG stream.
}

\usage{setSeedTP(x)}

\arguments{
	\item{x}{either a NULL object or an object of class \sQuote{TPmsmSeed}.}
}

\details{
An object of class \sQuote{TPmsmSeed} can be obtained by a previous
call to function \code{setSeedTP}, usually \code{setSeedTP(NULL)} or \code{setSeedTP()}.
The object can be saved and used as input on a later call to function \code{setSeedTP}
effectively restoring the seed of each individual RNG stream.
An object of class \sQuote{TPmsmSeed} can be manipulated or defined with arbitrary seeds,
however such procedure is not recommended. It is strongly recommended to input
objects of class \sQuote{TPmsmSeed} obtained from previous calls to function \code{setSeedTP}.
A seed of choice can be defined by calling function \code{\link{setPackageSeedTP}}.
}

\value{
Invisibly returns an object of class \sQuote{TPmsmSeed}.
\sQuote{TPmsmSeed} objects are implemented as a list of RNG stream seeds.
}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014) TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Karl A. T., Eubank R., Milovanovic J., Reiser M., Young D. (2014) Using RngStreams for parallel random number generation in C++ and R. \emph{Computational Statistics} \bold{29}(5), 1301--1320. \doi{10.1007/s00180-014-0492-3}

L'Ecuyer, P. (1999) Good parameters and implementations for combined multiple recursive random number generators. \emph{Operations Research} \bold{47}(1), 159–-164. \doi{10.1287/opre.47.1.159}

L’Ecuyer P., Simard R., Chen E. J., Kelton W. D. (2002) An object-oriented random-number package with many long streams and substreams. \emph{Operations Research} \bold{50}(6), 1073–-1075. \doi{10.1287/opre.50.6.1073.358}
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\note{
Unlike function \code{\link{setPackageSeedTP}} this function doesn't
recreate the RNG streams each time it is called.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{setPackageSeedTP}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Generate bivariate survival data
survTP0 <- dgpTP(n=100, corr=1, dist="weibull", dist.par=c(2, 7, 2, 7),
model.cens="exponential", cens.par = 6, state2.prob=0.6)

# Save seed
seed <- setSeedTP()

# Generate bivariate survival data
survTP1 <- dgpTP(n=100, corr=1, dist="weibull", dist.par=c(2, 7, 2, 7),
model.cens="exponential", cens.par = 6, state2.prob=0.6)

# The objects should be different
all.equal(survTP0, survTP1)

# Restore seed
setSeedTP(seed)

# Generate bivariate survival data
survTP2 <- dgpTP(n=100, corr=1, dist="weibull", dist.par=c(2, 7, 2, 7),
model.cens="exponential", cens.par = 6, state2.prob=0.6)

# Both objects were computed from the same seed and should be equal
all.equal(survTP1, survTP2)

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{environment}
\keyword{utilities}
