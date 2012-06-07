\name{dgpTP}
\alias{dgpTP}
\title{Generates bivariate survival data}
\description{Generates bivariate censored gap times from some known copula functions.}
\usage{dgpTP(n, corr, dist, dist.par, model.cens, cens.par,
state2.prob, to.data.frame = FALSE)}

\arguments{
\item{n}{Sample size.}
\item{corr}{Correlation parameter. Possible values for the bivariate exponential distribution are between 0 (for independency) and 1. Any value between 0 (not included) and 1 is accepted for the bivariate Weibull distribution.}
\item{dist}{Distribution. Possible bivariate distributions are "exponential" and "weibull".}
\item{dist.par}{Vector of parameters for the allowed distributions. Two (rate) parameters for the bivariate exponential distribution and four (2 location parameters and 2 scale parameters) for the bivariate Weibull distribution. See details below.}
\item{model.cens}{Model for censorship. Possible values are "uniform" and "exponential".}
\item{cens.par}{Parameter for the censorship distribution. For censure model equal to "exponential" the argument \code{cens.par} must be greater than 0. For censure model equal to "uniform" the argument must be greater or equal than 0.}
\item{state2.prob}{The proportion of individuals that enter state 2.}
\item{to.data.frame}{If \code{TRUE} returns a data.frame, if \code{FALSE} returns an object of class \code{survTP}. Defaults to \code{FALSE}.}
}

\details{
The bivariate exponential distribution, also known as Farlie-Gumbel-Morgenstern distribution is given by
\deqn{F(x,y)=F_1(x)F_2(y)[1+\delta((1-F_1(x))(1-F_2(y)))]}{F(x,y)=F_1(x)F_2(y)[1+d((1-F_1(x))(1-F_2(y)))}
for \eqn{x\ge0}{x\ge0} and \eqn{y\ge0}{y\ge0}. Where the marginal distribution functions \eqn{F_1}{F_1} and \eqn{F_2}{F_2} are exponential with rate parameters \eqn{\lambda_1}{l_1} and \eqn{\lambda_2}{l_2} and correlation parameter \eqn{\delta}{d}, \eqn{-1 \le \delta \le 1}{-1 \le \delta \le 1}. \cr

The bivariate Weibull distribution with two-parameter marginal distributions. It's survival function is given by
\deqn{S(x,y)=P(X>x,Y>y)=e^{-[(\frac{x}{\theta_1})^\frac{\beta_1}{\delta}+(\frac{y}{\theta_2})^\frac{\beta_2}{\delta}]^\delta}}{S(x,y)=P(X>x,Y>y)=exp^(-[(x/t_1)^\beta_1/d+(y/t_2)^b_2/d]^d)}

Where \eqn{0 < \delta \le 1}{0 < d \le 1} and each marginal distribution has shape parameter \eqn{\beta_i}{b_i} and a scale parameter \eqn{\theta_i}{t_i}, \eqn{i = 1, 2}.
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

	Mark E. Johnson (1987). "Multivariate Statistical Simulation", John Wiley and Sons.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{corrTP}}.
}

\examples{
#Example for the bivariate Exponential distribution
dgpTP(n=100, corr=1, dist="exponential", dist.par=c(1, 1),
model.cens="uniform", cens.par=3, state2.prob=0.5, to.data.frame=TRUE)
#Example for the bivariate Weibull distribution
dgpTP(n=100, corr=1, dist="weibull", dist.par=c(2, 7, 2, 7),
model.cens="exponential", cens.par = 0.08, state2.prob=0.6)
}
