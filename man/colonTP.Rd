\encoding{UTF-8}
\name{colonTP}
\alias{colonTP}
\docType{data}
\title{Chemotherapy for Stage B/C colon cancer}
\description{These are data from one of the first successful trials of
	adjuvant chemotherapy for colon cancer. Levamisole is a low-toxicity
	compound previously used to treat worm infestations in animals; 5-FU
	is a moderately toxic (as these things go) chemotherapy agent.}
\usage{data(colonTP)}

\format{
	A data frame with 929 observations on the following 15 variables.
	\describe{
		\item{\code{time1}}{Time to recurrence/censoring/death, whichever occurs first.}
		\item{\code{event1}}{Recurrence/censoring/death indicator (recurrence/dead=1, alive=0).}
		\item{\code{Stime}}{Time to death/censoring.}
		\item{\code{event}}{Censoring indicator (dead=1, alive=0).}
		\item{\code{rx}}{Treatment - Obs(ervation), Lev(amisole), Lev(amisole)+5-FU.}
		\item{\code{sex}}{1=male.}
		\item{\code{age}}{Age in years.}
		\item{\code{obstruct}}{Obstruction of colon by tumour.}
		\item{\code{perfor}}{Perforation of colon.}
		\item{\code{adhere}}{Adherence to nearby organs.}
		\item{\code{nodes}}{Number of lymph nodes with detectable cancer.}
		\item{\code{differ}}{Differentiation of tumour (1=well, 2=moderate, 3=poor).}
		\item{\code{extent}}{Extent of local spread (1=submucosa, 2=muscle, 3=serosa, 4=contiguous structures).}
		\item{\code{surg}}{Time from surgery to registration (0=short, 1=long).}
		\item{\code{node4}}{More than 4 positive lymph nodes.}
	}
}

\note{The study is originally described in Laurie (1989).
	The main report is found in Moertel (1990).  This data set is closest
	to that of the final report in Moertel (1991).
	A version of the data with less follow-up time was used in the
	paper by Lin (1994).
}

\references{
	JA Laurie, CG Moertel, TR Fleming, HS Wieand, JE Leigh, J Rubin,
	GW McCormack, JB Gerstner, JE Krook and J Malliard.  Surgical
	adjuvant therapy of large-bowel carcinoma: An evaluation of
	levamisole and the combination of levamisole and fluorouracil:
	The North Central Cancer Treatment Group and the Mayo Clinic.
	\emph{Journal of Clinical Oncology}, 7:1447-1456, 1989.

	DY Lin.  Cox regression analysis of multivariate failure time data:
	the marginal approach.  \emph{Statistics in Medicine}, 13:2233-2247, 1994.

	CG Moertel, TR Fleming, JS MacDonald,
	DG Haller, JA Laurie, PJ Goodman, JS Ungerleider,
	WA Emerson, DC Tormey, JH Glick, MH Veeder and JA Maillard.
	Levamisole and fluorouracil for adjuvant therapy of
	resected colon carcinoma. \emph{New England Journal of Medicine}, 332:352-358,
	1990.

	CG Moertel,  TR Fleming,  JS MacDonald, DG Haller, JA Laurie, CM Tangen,
	JS Ungerleider, WA Emerson, DC Tormey, JH Glick, MH Veeder and JA
	Maillard, Fluorouracil plus Levamisole as and effective adjuvant
	therapy after resection of stage II colon carcinoma: a final report.
	\emph{Annals of Internal Medicine}, 122:321-326, 1991.
}

\keyword{datasets}
