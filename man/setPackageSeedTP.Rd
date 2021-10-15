\encoding{UTF-8}
\name{setPackageSeedTP}
\alias{setPackageSeedTP}
\title{Set the initial package seed}

\description{
The random number generator (RNG) with multiple independent streams developed by
L'Ecuyer et al. (2002) is used for parallel computation of uniform pseudorandom numbers.
Package \pkg{TPmsm} makes extensive use of uniform pseudorandom numbers,
particularly for the bootstrapping statistical techniques and for
the generation of univariate and multivariate pseudorandom data.
This function defines the seed for the creation of RNG streams.
}

\usage{setPackageSeedTP(seed=12345)}

\arguments{
	\item{seed}{A vector of one to six integers.
		Defaults to \code{rep(x=12345, times=6)}.
	}
}

\details{
If the user defines a vector with length lower than six
as seed, then the seed is internally defined as a vector
of length six with the first elements equal to the user
defined values, and the leaving elements equal to \code{12345}.
If a vector with more than six elements is provided as seed,
then only the first six elements are used.
}

\value{Invisibly returns NULL.}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Karl A. T., Eubank R., Milovanovic J., Reiser M., Young D. (2014). Using RngStreams for parallel random number generation in C++ and R. \emph{Computational Statistics}, \bold{29}(5), 1301-1320. \doi{10.1007/s00180-014-0492-3}

L'Ecuyer, P. (1999). Good parameters and implementations for combined multiple recursive random number generators. \emph{Operations Research}, \bold{47}(1), 159–-164. \doi{10.1287/opre.47.1.159}

L’Ecuyer P., Simard R., Chen E. J., Kelton W. D. (2002). An object-oriented random-number package with many long streams and substreams. \emph{Operations Research}, \bold{50}(6), 1073–-1075. \doi{10.1287/opre.50.6.1073.358}
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\note{
When package \pkg{TPmsm} loads, an initial set of RNG streams
is created, one stream for each thread available for parallel
computation. The initial set of RNG streams is created from the
package seed \code{c(12345, 12345, 12345, 12345, 12345, 12345)}.
Every time this function is called, the old set of RNG streams
is deleted, and a new set of RNG streams is created from the
user defined package seed. After the creation of each new RNG stream,
the internally stored package seed changes. So each RNG stream is
created from a different package seed, and yields different sets of
pseudorandom numbers.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{setSeedTP}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2)

# Define package seed
seed <- rep(x=1, times=6)

# Set package seed
setPackageSeedTP(seed)

# Create survTP object
data(heartTP)
heartTP_obj <- with(heartTP, survTP(time1, event1, Stime, event))

# Compute transition probabilities with confidence band
TPmsm0 <- transAJ(object=heartTP_obj, s=33, t=412, conf=TRUE,
conf.level=0.9, method.boot="percentile")

# Compute transition probabilities with confidence band
TPmsm1 <- transAJ(object=heartTP_obj, s=33, t=412, conf=TRUE,
conf.level=0.9, method.boot="percentile")

# The objects should be different
all.equal(TPmsm0, TPmsm1)

# Set package seed
setPackageSeedTP(seed)

# Compute transition probabilities with confidence band
TPmsm2 <- transAJ(object=heartTP_obj, s=33, t=412, conf=TRUE,
conf.level=0.9, method.boot="percentile")

# Both objects were computed from the same seed and should be equal
all.equal(TPmsm0, TPmsm2)

# Restore the number of threads
setThreadsTP(nth)
}

\keyword{environment}
\keyword{utilities}
