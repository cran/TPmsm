\encoding{UTF-8}
\name{TPmsm-package}
\alias{TPmsm-package}
\alias{TPmsm}
\docType{package}
\title{Estimation of transitions probabilities in multistate models}
\description{
	The \pkg{TPmsm} software contains functions that compute estimates for the transition probabilities in the
	illness-death model and or the three-state progressive model. This package implements seven different estimators.
	Being five of them non-parametric and two of them semi-parametric (PAJ and KMPW). The implemented estimators are the
	Aalen-Johansen estimator (AJ), Presmoothed Aalen-Johansen estimator (PAJ), Kaplan-Meier Weighted estimator (KMW),
	Presmoothed Kalpan-Meier Weighted estimator (KMPW), Inverse Probability of Censoring estimator (IPCW),
	Lin estimator (LIN) and Location-Scale estimator (LS). The Inverse Probability of Censoring (IPCW) and Lin (LIN)
	estimators also permit to compute transition probabilities conditioned on a single covariate. Bootstrap confidence
	bands can be computed for each of the mentioned estimators. Several graphical plots of the transition probabilities
	with or without confidence bands can be drawn. To aid in the study of the statistical properties of the implemented
	estimators, functions to generate pseudo-random data for some well-known multivariate distributions were implemented.
}

\details{
	\tabular{ll}{
		Package: \tab TPmsm\cr
		Type: \tab Package\cr
		Version: \tab 1.1.1\cr
		Date: \tab 2014-02-06\cr
		License: \tab GPL (>= 2)\cr
		LazyLoad: \tab yes\cr
		LazyData: \tab yes\cr
	}
}

\author{
	Artur Agostinho Araújo, Javier Roca-Pardiñas \email{roca@uvigo.es}\cr
	and Luís Meira Machado \email{lmachado@math.uminho.pt}\cr
	Maintainer: Artur Agostinho Araújo \email{artur.stat@gmail.com}
}

\references{
	Aalen O. O., Johansen S. (1978) An Empirical Transition Matrix for Nonhomogeneous Markov Chains Based on Censored Observations. \emph{Scandinavian Journal of Statistics}, \bold{5}, 141-150.

	Allignol A., Schumacher M., Beyersmann J. (2011). Empirical Transition Matrix of Multi-State Models: The etm Package. \emph{Journal of Statistical Software}, \bold{38(4)}, 1-15.

	Amorim A. P., de Uña-Álvarez J., Meira Machado L. F. (2011). Presmoothing the transition probabilities in the illness-death model, \emph{Statistics and Probability Letters}, \bold{81}, 797-806.

	Davison, A. C. and Hinkley, D. V. (1997) Bootstrap Methods and their Application, Chapter 5. Cambridge University Press.

	Johnson M. E. (1987). Multivariate Statistical Simulation, John Wiley and Sons.

	Johnson N., Kotz S. (1972). Distributions in statistics: continuous multivariate distributions, John Wiley and Sons.

	Lu J., Bhattacharya G. (1990). Some new constructions of bivariate weibull models, \emph{Annals of Institute of Statistical Mathematics}, \bold{42}, 543-559.

	Meira Machado L. F., de Uña-Álvarez J., Cadarso-Suárez C. (2006). Nonparametric estimation of transition probabilities in a non-Markov illness-death model, \emph{Lifetime Data Anal}, \bold{12}, 325-344.

	Meira-Machado L., de Uña-Álvarez J. and Datta S. Conditional Transition Probabilities in a non-Markov Illness-death Model. Discussion Papers in Statistics and Operations Research n 11/03, 2011. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402 - 2007). This file can be downloaded from: \url{http://webs.uvigo.es/depc05/reports/12_05.pdf}

	Meira-Machado L., Roca-Pardiñas J. (2011). p3state.msm: Analyzing Survival Data from an Illness-Death Model. \emph{Journal of Statistical Software}, \bold{38(3)}, 1-18.

	Meira-Machado L., Roca-Pardiñas J., Van Keilegom I., and Cadarso-Suárez C. (2013). Bandwidth Selection for the Estimation of Transition Probabilities in the Location-Scale Progressive Three-State Model. \emph{Computational Statistics}, \bold{28}, 2185-2210.

	Meira-Machado L., Roca-Pardiñas J., Van Keilegom I. and Cadarso-Suárez C. Estimation of transition probabilities in a non-Markov model with successive survival times. Discussion paper 2010. This file can be downloaded from: \url{http://sites.uclouvain.be/IAP-Stat-Phase-V-VI/ISBApub/dp2010/DP1053.pdf}

	Moreira A., de Uña-Álvarez J. and Meira-Machado L. Presmoothing the Aalen-Johansen estimator of transition probabilities. Discussion Papers in Statistics and Operation Research n 11/03, 2011. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402 - 2007). This file can be downloaded from: \url{http://webs.uvigo.es/depc05/reports/11_03.pdf}

	Van Keilegom I., de Uña-Álvarez J. and Meira-Machado L. (2011) Nonparametric location-scale models for successive survival times under dependent censoring. \emph{Journal of Statistical Planning and Inference}, \bold{141}, 1118-1131.
}

\keyword{package}
