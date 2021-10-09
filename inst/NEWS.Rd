\name{NEWS}

\title{NEWS file for the TPmsm package}

\section{Changes in version 1.2.4}{
	\itemize{
		\item{Implemented vignette compression.}
		\item{Improved bibliographic references.}
	}
}

\section{Changes in version 1.2.3}{
	\itemize{
		\item{Added package vignette.}
		\item{Improved the documentation.}
	}
}

\section{Changes in version 1.2.2}{
	\itemize{
		\item{Fixed compiler warning in "RngStream.c" source file.}
		\item{Temporarily disabled OpenMP with Clang compiler in
		"transLS.c" source file.}
	}
}

\section{Changes in version 1.2.1}{
	\itemize{
		\item{Corrected DESCRIPTION and NAMESPACE files.}
	}
}

\section{Changes in version 1.2.0}{
	\itemize{
		\item{Improved the documentation.}
		\item{Renamed threadsTP function to setThreadsTP.}
		\item{Implemented RngStreams C package for
		parallel random number generation.}
		\item{Added setPackageSeedTP function for
		defining the seed of the RngStreams C package.}
		\item{Added setSeedTP function for saving and restoring
		the individual seeds of the RNG streams used by the
		RngStreams C package.}
		\item{Added research article to the CITATION file.}
	}
}

\section{Changes in version 1.1.1}{
	\itemize{
		\item{Improved the documentation.}
		\item{Removed unnecessary require calls from R code.}
	}
}

\section{Changes in version 1.1.0}{
	\itemize{
		\item{Changed package encoding to UTF-8.}
		\item{Added lines methods for TPmsm and TPCmsm objects.}
		\item{Added CITATION file.}
		\item{Changed parameterization of the bivariate and univariate
		exponential models of distribution available for the dgpTP function.
		Now the mean of the univariate exponential model of distribution is
		equal to it's parameter.}
		\item{Removed symmetric and asymmetric options from window
		parameter in transLS, transLIN and transIPCW functions.}
		\item{Changes were made to the kernel algorithms that might
		improve speed up to around 25\% in bootstrap. This is particularly
		evident in the cross-validation algorithm available in the transLS
		function.}
		\item{Hided unneeded foreign language functions inside the shared
		object or dynamic link library.}
		\item{Added option to compute the bandwidth by cross-validation
		for both the location and scale functions. This applies to
		the Location-Scale estimator.}
		\item{Fixed a bug in both the AJ and PAJ algorithms that might result
		in wrong transition probability estimates from state 1 to state 2
		and from state 1 to state 3.}
		\item{Added as.data.frame.survTP function.}
		\item{Renamed adaptp3state function to TPmsmOut. Added option
		to convert a data.frame in the "TPmsm" format to a data.frame
		in the "etm" format.}
		\item{Removed to.data.frame argument from dgpTP function.
		Now dgpTP always returns an object of class survTP.}
		\item{Transition probability bootstrap confidence bands from
		state 1 to state 2 and from state 1 to state 3 could be wrong
		in previous versions. This situation affected all estimators.
		Fixed in this version.}
		\item{Improved the documentation.}
		\item{Added threadsTP function.}
		\item{Changed some examples to allow them to run faster.
		Now all examples run with a maximum of two threads.}
	}
}

\section{Changes in version 1.0.1}{
	\itemize{
		\item{Fixed failed compilation under older releases of R.}
	}
}

\section{Changes in version 1.0}{
	\itemize{
		\item{First release.}
	}
}
