\name{selectionMet}
\alias{selectionMet}
\title{
Metabolites selection
}
\description{
This function select metabolites from the list generated by \code{\link{readMet}} function.
}
\usage{
selectionMet(doc, sel)
}
\arguments{
  \item{doc}{A list of metabolites information produced by \code{\link{readMet}} function.}
  \item{sel}{A vector of metabolite's HMDBcode that will be selected}
}

\value{
 \item{doc}{A doc list contains only the selcted metabolites.}
}

\seealso{
  \code{\link{readMet}}, \code{\link{nameMet}}
}


\examples{

data(Metabolites)
doc=selectionMet(Metabolites$readMet,c("HMDB0000299","HMDB0000881"))
nameMet(doc)
}

