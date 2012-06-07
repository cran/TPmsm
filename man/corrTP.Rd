\name{corrTP}
\alias{corrTP}
\title{Correlation between two gap times}
\description{Provides the correlation between the bivariate times for some copula distributions.}
\usage{corrTP(dist, corr, dist.par)}

\arguments{
\item{dist}{The distribution. Possible bivariate distributions are "exponential" and "weibull".}
\item{corr}{Correlation parameter. Possible values for the bivariate exponential distribution are between -1 and 1 (0 for independency). Any value between 0 (not included) and 1 is accepted for the bivariate Weibull distribution.}
\item{dist.par}{Vector of parameters for the allowed distributions. Two (rate) parameters for the bivariate exponential distribution and four (2 location parameters and 2 scale parameters) for the bivariate Weibull distribution. See details below.}
}

\details{
The bivariate exponential distribution, also known as Farlie-Gumbel-Morgenstern distribution is given by
\deqn{F(x,y)=F_1(x)F_2(y)[1+\delta((1-F_1(x))(1-F_2(y)))]}{F(x,y)=F_1(x)F_2(y)[1+\delta((1-F_1(x))(1-F_2(y)))}
for \eqn{x\ge0}{x\ge0} and \eqn{y\ge0}{y\ge0}. Where the marginal distribution functions \eqn{F_1}{F_1} and \eqn{F_2}{F_2} are exponential with rate parameters \eqn{\lambda_1}{\lambda_1} and \eqn{\lambda_2}{\lambda_2} and correlation parameter \eqn{\delta}{\delta}, \eqn{-1 \le \delta \le 1}{-1 \le \delta \le 1}. \cr

The bivariate Weibull distribution with two-parameter marginal distributions. It's survival function is given by
\deqn{S(x,y)=P(X>x,Y>y)=e^{-[(\frac{x}{\theta_1})^\frac{\beta_1}{\delta}+(\frac{y}{\theta_2})^\frac{\beta_2}{\delta}]^\delta}}{S(x,y)=P(X>x,Y>y)=exp^(-[(x/\theta_1)^\beta_1/\delta+(y/\theta_2)^\beta_2/\delta]^\delta)}

Where \eqn{0 < \delta \le 1}{0 < \delta \le 1} and each marginal distribution has shape parameter \eqn{\beta_i}{\beta_i} and a scale parameter \eqn{\theta_i}{\theta_i}, \eqn{i = 1, 2}.
}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

\author{Artur Agostinho Araujo, Javier Roca-Pardinas and Luis Meira Machado}

\references{
	Johnson N, Kotz S (1972). "Distribution in statistics: continuous multivariate distributions", John Wiley and Sons.

	Lu J, Bhattacharya G (1990). "Some new constructions of bivariate weibull models", Annals of Institute of Statistical Mathematics, \bold{42}, 543-559.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dgpTP}}.
}

\examples{
#Example for the bivariate Weibull distribution
corrTP(dist = "weibull", corr = 0.5, dist.par = c(2, 7, 2, 7))
#Example for the bivariate Exponential distribution
corrTP(dist = "exponential", corr = 1, dist.par = c(1, 1))
}
