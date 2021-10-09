\encoding{UTF-8}
\name{colonTP}
\alias{colonTP}
\docType{data}
\title{Chemotherapy for Stage B/C colon cancer}

\description{
These are data from one of the first successful trials of
adjuvant chemotherapy for colon cancer. Levamisole is a low-toxicity
compound previously used to treat worm infestations in animals; 5-FU
is a moderately toxic (as these things go) chemotherapy agent.
}

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

\note{
The study is originally described in Laurie (1989).
The main report is found in Moertel (1990).  This data set is closest
to that of the final report in Moertel (1991).
A version of the data with less follow-up time was used in the
paper by Lin (1994).
}

\references{
Laurie J. A., Moertel C. G., Fleming T. R., Wieand H. S., Leigh J. E., Rubin J., McCormack G. W., Gerstner J. B., Krook J. E., Malliard J. (1989) Surgical adjuvant therapy of large-bowel carcinoma: An evaluation of levamisole and the combination of levamisole and fluorouracil: The North Central Cancer Treatment Group and the Mayo Clinic. \emph{Journal of Clinical Oncology} \bold{7}(10), 1447--1456. \doi{10.1200/JCO.1989.7.10.1447}

Lin D. Y. (1994) Cox regression analysis of multivariate failure time data: the marginal approach. \emph{Statistics in Medicine} \bold{13}(21), 2233--2247. \doi{10.1002/sim.4780132105}

Moertel C. G., Fleming T. R., MacDonald J. S., Haller D. G., Laurie J. A., Goodman P.J., Ungerleider J.S., Emerson W.A., Tormey D.C., Glick J.H., Veeder M.H., Maillard J.A. (1990) Levamisole and fluorouracil for adjuvant therapy of resected colon carcinoma. \emph{New England Journal of Medicine} \bold{322}(6), 352--358. \doi{10.1056/NEJM199002083220602}

Moertel C. G., Fleming T. R., MacDonald J. S., Haller D. G., Laurie J. A., Tangen C. M., Ungerleider J. S., Emerson W. A., Tormey D. C., Glick J. H., Veeder M. H., Maillard J. A. (1995) Fluorouracil plus Levamisole as an effective adjuvant therapy after resection of stage II colon carcinoma: a final report. \emph{Annals of Internal Medicine} \bold{122}(5), 321--326. \doi{10.7326/0003-4819-122-5-199503010-00001}
}

\keyword{datasets}
