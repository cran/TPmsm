\encoding{UTF-8}
\name{contour.TPCmsm}
\alias{contour.TPCmsm}
\title{contour method for a TPCmsm object}

\description{
contour method for an object of class \sQuote{TPCmsm}.
Creates a contour plot of the transition probabilities.
}

\usage{\S3method{contour}{TPCmsm}(x, contour.type="tc", tr.choice, nlevels=20, levels=pretty(zlim, nlevels),
xlim, ylim, zlim=c(0, 1), col=grey(0.4), xlab, ylab, main="", sub="",
add=FALSE, las=1, conf.int=FALSE, legend=TRUE, curvlab, ...)}

\arguments{
	\item{x}{An object of class \sQuote{TPCmsm}.}
	\item{contour.type}{A character string specifying the type of contour.
		If \dQuote{tc} the contour with time in the x axis, covariate in the y axis and transition probability in the z axis is drawn.
		If \dQuote{ct} the contour with covariate in the x axis, time in the y axis and transition probability in the z axis is drawn.
		Defaults to \dQuote{tc}.
	}
	\item{tr.choice}{Character vector of the form \sQuote{c(\dQuote{from to}, \dQuote{from to})}
		specifying which transitions should be plotted.
		Default, all the transition probabilities are plotted.
	}
	\item{nlevels}{The number of levels to divide the range of z. Defaults to 20 levels.}
	\item{levels}{Numeric vector of levels at which to draw contour lines.
		Defaults to \code{pretty(zlim, nlevels)}.
	}
	\item{xlim}{Limits of x-axis for the plot.}
	\item{ylim}{Limits of y-axis for the plot.}
	\item{zlim}{Limits of z-axis for the plot. Defaults to \code{c(0, 1)}.}
	\item{col}{Color for the lines drawn. Defaults to \code{grey(0.4)}.}
	\item{xlab}{x-axis label. If \code{contour.type="tc"} defaults to \dQuote{Time}.
		If \code{contour.type="ct"} defaults to \dQuote{Covariate}.
	}
	\item{ylab}{y-axis label. If \code{contour.type="tc"} defaults to \dQuote{Covariate}.
		If \code{contour.type="ct"} defaults to \dQuote{Time}.
	}
	\item{main}{The main title for the plot. By default no main title is added.}
	\item{sub}{A sub title for the plot. By default no sub title is added.}
	\item{add}{logical. If TRUE, add to a current plot.}
	\item{las}{The style of labeling to be used. The default is to use horizontal labeling.}
	\item{conf.int}{Logical. Whether to display contour plots of confidence regions. Default is FALSE.}
	\item{legend}{A logical specifying if a legend should be added.}
	\item{curvlab}{A character or expression vector to appear in the legend.
		Default is the name of the transitions.
	}
	\item{\dots}{Further arguments for contour.}
}

\value{
No value is returned.
}

\author{Artur Araújo, Javier Roca-Pardiñas and Luís Meira-Machado}

\note{
The device is divided by the number of transitions specified by argument \code{tr.choice}.
Being the number of columns equal to the number of transitions.
If argument \code{conf.int=TRUE} the device is further divided to make room for the confidence regions.
In this case two rows are added, one for each side of the confidence region.
So if \code{conf.int=TRUE} the center row provides the contour of the estimates.
The upper row provides the upper side of the confidence region.
And the lower row provides the lower side of the confidence region.
}

\references{
Araújo A, Meira-Machado L, Roca-Pardiñas J (2014). TPmsm: Estimation of the Transition Probabilities in
3-State Models. \emph{Journal of Statistical Software}, \bold{62}(4), 1-29. \doi{10.18637/jss.v062.i04}

Meira-Machado L., de Uña-Álvarez J., Datta S. (2011). Conditional Transition Probabilities in a non-Markov Illness-death Model. Discussion Papers in Statistics and Operations Research n 11/03. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402-2007). \url{https://depc05.webs.uvigo.es/reports/12_05.pdf}
}

\seealso{
\code{\link{contour}}.
}

\examples{
# Set the number of threads
nth <- setThreadsTP(2);

# Create survTP object
data(colonTP);
colonTP_obj <- with( colonTP, survTP(time1, event1, Stime, event, age=age) );

# Compute IPCW conditional transition probabilities without confidence band
TPCmsm_obj <- transIPCW(colonTP_obj, s=57, t=310, x=0);

# Plot contour with Time in the x-axis
contour(TPCmsm_obj, contour.type="tc", tr.choice=c("1 1", "1 2", "2 2"), ylab="Age");

# Plot contour with Time in the y-axis
contour(TPCmsm_obj, contour.type="ct", tr.choice=c("1 1", "1 2", "1 3"), xlab="Age");

# Restore the number of threads
setThreadsTP(nth);
}

\keyword{aplot}
\keyword{hplot}
\keyword{methods}
\keyword{multivariate}
\keyword{survival}
