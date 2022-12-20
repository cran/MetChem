\docType{data}

\name{Metabolites}

\alias{Metabolites}

\title{Metabolomic Dataset}

\description{
This dataset consists of a list of the metabolites as returned by the function \code{readMet} and concentration value of each metabolites.
}

\usage{data(Metabolites)}

\value{
A list with the following elements in the variable \code{Metabolites}:
	\item{concentration}{A matrix containing the concentration of each metabolites.}
	\item{name}{A vector of the HMDB IDs of the each metabolites .}
	\item{SMILE}{A vector of secondary IDs of the each metabolites .}
	\item{HMDB}{A vector containing all information relative to the feature of the each metabolite.}
	\item{ReadMet}{A list of metabolites information produced by \code{\link{readMet}} function.}
	}


\examples{
 data(Metabolites)

}

\keyword{datasets}
