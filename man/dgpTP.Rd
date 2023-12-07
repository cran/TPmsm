\encoding{UTF-8}
\name{dgpTP}
\alias{dgpTP}
\title{Generates bivariate survival data}

\description{Generates bivariate censored gap times from some known copula functions.}

\usage{dgpTP(n, corr, dist, dist.par, model.cens, cens.par, state2.prob)}

\arguments{
	\item{n}{Sample size.}
	\item{corr}{Correlation parameter.
		Possible values for the bivariate exponential distribution are between -1 and 1 (0 for independency).
		Any value between 0 (not included) and 1 (1 for independency) is accepted for the bivariate Weibull distribution.
	}
	\item{dist}{Distribution. Possible bivariate distributions are \dQuote{exponential} and \dQuote{weibull}.}
	\item{dist.par}{Vector of parameters for the allowed distributions.
		Two (scale) parameters for the bivariate exponential distribution and four (2 location parameters and 2 scale parameters)
		for the bivariate Weibull distribution. See details below.
	}
	\item{model.cens}{Model for censorship. Possible values are \dQuote{uniform} and \dQuote{exponential}.}
	\item{cens.par}{Parameter for the censorship distribution.
		For censure model equal to \dQuote{exponential} the argument \code{cens.par} must be greater than 0.
		For censure model equal to \dQuote{uniform} the argument must be greater or equal than 0.
	}
	\item{state2.prob}{The proportion of individuals that enter state 2.}
}

\details{
The bivariate exponential distribution, also known as Farlie-Gumbel-Morgenstern distribution is given by
\deqn{F(x,y)=F_1(x)F_2(y)[1+\alpha(1-F_1(x))(1-F_2(y))]}{F(x,y)=F_1(x)F_2(y)[1+\alpha(1-F_1(x))(1-F_2(y))]}
for \eqn{x\ge0}{x\ge0} and \eqn{y\ge0}{y\ge0}. Where the marginal distribution functions \eqn{F_1}{F_1} and \eqn{F_2}{F_2} are exponential with scale parameters \eqn{\theta_1}{\theta_1} and \eqn{\theta_2}{\theta_2} and correlation parameter \eqn{\alpha}{\alpha}, \eqn{-1 \le \alpha \le 1}{-1 \le \alpha \le 1}. \cr

The bivariate Weibull distribution with two-parameter marginal distributions. It's survival function is given by
\deqn{S(x,y)=P(X>x,Y>y)=e^{-[(\frac{x}{\theta_1})^\frac{\beta_1}{\delta}+(\frac{y}{\theta_2})^\frac{\beta_2}{\delta}]^\delta}}{S(x,y)=P(X>x,Y>y)=exp^(-[(x/\theta_1)^(\beta_1/\delta)+(y/\theta_2)^(\beta_2/\delta)]^\delta)}

Where \eqn{0 < \delta \le 1}{0 < \delta \le 1} and each marginal distribution has shape parameter \eqn{\beta_i}{\beta_i} and a scale parameter \eqn{\theta_i}{\theta_i}, \eqn{i = 1, 2}.
}

\value{An object of class \sQuote{survTP}.}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Devroye L. (1986). \emph{Non-Uniform Random Variate Generation}, New York: Springer-Verlag.

Johnson M. E. (1987). \emph{Multivariate Statistical Simulation}, John Wiley and Sons.

Johnson N., Kotz S. (1972). \emph{Distributions in statistics: continuous multivariate distributions}, John Wiley and Sons.

Lu J., Bhattacharya G. (1990). Some new constructions of bivariate weibull models. \emph{Annals of Institute of Statistical Mathematics}, \bold{42}(3), 543-559. \doi{10.1007/BF00049307}

Meira-Machado L., Faria S. (2014). A simulation study comparing modeling approaches in an illness-death multi-state model. \emph{Communications in Statistics - Simulation and Computation}, \bold{43}(5), 929-946. \doi{10.1080/03610918.2012.718841}

Meira-Machado, L., Sestelo M. (2019). Estimation in the progressive illness-death model: a nonexhaustive
review. \emph{Biometrical Journal}, \bold{61}(2), 245–263. \doi{10.1002/bimj.201700200}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{corrTP}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2);

# Example for the bivariate Exponential distribution
dgpTP(n=100, corr=1, dist="exponential", dist.par=c(1, 1),
model.cens="uniform", cens.par=3, state2.prob=0.5);

# Example for the bivariate Weibull distribution
dgpTP(n=100, corr=1, dist="weibull", dist.par=c(2, 7, 2, 7),
model.cens="exponential", cens.par = 6, state2.prob=0.6);

# Restore the number of threads
setThreadsTP(nth);
}

\keyword{datagen}
\keyword{distribution}
\keyword{multivariate}
\keyword{survival}
