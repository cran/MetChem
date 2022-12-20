\name{features}
\alias{features}
\title{
Cluster features extraction 
}
\description{
This function finds features associated with each cluster.
}
\usage{

features(doc,cla,clustering)
}
\arguments{
  \item{doc}{The output of the \code{\link{readMet}} function.}
  \item{cla}{The output of \code{\link{diseasesMet}}, \code{\link{enzymesMet}}, \code{\link{pathwaysMet}}, \code{\link{propertiesMet}}, \code{\link{substituentsMet}}, or \code{\link{taxonomyMet}}  functions.}
  \item{clustering}{A clustering vector of the metabolites.}
}

\value{
  A list of p-value calculated using Fisher test for cluster associted features.
}

\seealso{
  \code{\link{KODAMA.chem.sim}}, \code{\link{tree.cutting}}, \code{\link{substituentsMet}}
}

\examples{
\donttest{
data(Metabolites)
res=KODAMA.chem.sim(Metabolites$SMILE)
clu=tree.cutting(res,max_nc = 30)
cla=substituentsMet(Metabolites$readMet)
f=features(Metabolites$readMet,cla,clu$clusters[,"Clusters 23"])
}
}
