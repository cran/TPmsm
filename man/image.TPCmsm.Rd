\encoding{UTF-8}
\name{image.TPCmsm}
\alias{image.TPCmsm}
\title{image method for a TPCmsm object}
\description{
	image method for an object of class 'TPCmsm'.
	Creates a grid of colored or gray-scale rectangles with colors corresponding to the values of the transition probabilities.
}
\usage{\S3method{image}{TPCmsm}(x, image.type="tc", tr.choice, xlim, ylim, zlim=c(0, 1), col, xlab, ylab,
main, sub, key.title, key.axes, las=1, conf.int=FALSE, legend=TRUE, curvlab,
contour=TRUE, nlevels=20, levels=pretty(zlim, nlevels), ...)}

\arguments{
	\item{x}{An object of class 'TPCmsm'.}
	\item{image.type}{A character string specifying the type of image.
		If 'tc' the image with time in the x axis, covariate in the y axis and transition probability in the z axis is drawn.
		If 'ct' the image with covariate in the x axis, time in the y axis and transition probability in the z axis is drawn.
		Defaults to 'tc'.}
	\item{tr.choice}{Character vector of the form 'c("from to", "from to")' specifying which transitions should be plotted.
		Default, all the transition probabilities are plotted.}
	\item{xlim}{Limits of x-axis for the plot.}
	\item{ylim}{Limits of y-axis for the plot.}
	\item{zlim}{Limits of z-axis for the plot. Defaults to \code{c(0, 1)}.}
	\item{col}{Vector of colour. Defaults to \code{heat.colors(nlevels)[nlevels:1]}.}
	\item{xlab}{x-axis label. If \code{image.type="tc"} defaults to "Time". If \code{image.type="ct"} defaults to "Covariate".}
	\item{ylab}{y-axis label. If \code{image.type="tc"} defaults to "Covariate". If \code{image.type="ct"} defaults to "Time".}
	\item{main}{The main title for the plot. By default no main title is added.}
	\item{sub}{A sub title for the plot. By default no sub title is added.}
	\item{key.title}{Statements which add titles for the plot key.}
	\item{key.axes}{Statements which draw axes on the plot key. This overrides the default axis.}
	\item{las}{The style of labeling to be used. The default is to use horizontal labeling.}
	\item{conf.int}{Logical. Whether to display images of confidence regions. Default is FALSE.}
	\item{legend}{A logical specifying if a legend should be added.}
	\item{curvlab}{A character or expression vector to appear in the legend. Default is the name of the transitions.}
	\item{contour}{If \code{TRUE} contour lines are added to the image. Default is \code{TRUE}.}
	\item{nlevels}{The number of levels to divide the range of z. Defaults to 20 levels.}
	\item{levels}{Numeric vector of levels at which to draw contour lines. Defaults to
		\code{pretty(zlim, nlevels)}. The breaks of image are set equal to levels.}
	\item{\dots}{Further arguments for image.}
}

\value{
	No value is returned.
}

\author{Artur Agostinho Araújo, Javier Roca-Pardiñas and Luís Meira Machado}

\note{
	The device is divided by the number of transitions specified by argument \code{tr.choice}.
	Being the number of columns equal to the number of transitions.
	If argument \code{conf.int=TRUE} the device is further divided to make room for the confidence regions.
	In this case two rows are added, one for each side of the confidence region.
	So if \code{conf.int=TRUE} the center row provides the image of the estimates.
	The upper row provides the upper side of the confidence region.
	And the lower row provides the lower side of the confidence region.
}

\references{
	Meira-Machado L., de Uña-Álvarez J. and Datta S. Conditional Transition Probabilities in a non-Markov Illness-death Model. Discussion Papers in Statistics and Operation Research n 11/03, 2011. Department of Statistics and Operations Research, University of Vigo (ISSN: 1888-5756, Deposito Legal VG 1402 - 2007). This file can be downloaded from: \url{http://webs.uvigo.es/depc05/reports/12_05.pdf}
}

\seealso{
	\code{\link{contour}},
	\code{\link{image}}.
}

\examples{
# set the number of threads
nth <- threadsTP(2)

# Create survTP object
data(heartTP)
heartTP_obj <- with( heartTP, survTP(time1, event1, Stime, event, age=age) )

# Compute LIN conditional transition probabilities with confidence band
TPCmsm_obj <- transLIN(heartTP_obj, s=57, t=310, x=0, conf=TRUE, n.boot=100,
method.boot="basic")

# Plot image with Time in the x-axis
image(TPCmsm_obj, image.type="tc", tr.choice=c("1 1", "1 2", "2 2"), conf.int=TRUE,
ylab="Age")

# Plot image with Time in the y-axis
image(TPCmsm_obj, image.type="ct", tr.choice=c("1 1", "1 2", "1 3"), conf.int=TRUE,
xlab="Age")

# restore the number of threads
threadsTP(nth)
}

\keyword{hplot}
\keyword{methods}
\keyword{multivariate}
\keyword{survival}
