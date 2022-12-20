\name{WMCSA}
\alias{WMCSA}
\title{
Weighted Metabolite Chemical Structural Analysis
}
\description{
Summarize metabolites concetration in each of identified clusters using the module eigenvalue (eigen-metabolite) for calculating module membership measures. 
}
\usage{
WMCSA(data,clu,nclusters)
}

\arguments{
\item{data}{
dataset of different metabolite concentration in differnt samples.
  }
\item{clu}{
a list in the same format of the output of the function \code{\link{tree.cutting}} function. 
}
\item{nclusters}{
number of clusters.
}
}

\value{
This function returns a matrix as output represent similarity score of metabolites within the same cluster among different samples.

}

\seealso{\code{\link{KODAMA.chem.sim}}, \code{\link{tree.cutting}}
}
\examples{
\donttest{

data(Metabolites)

res=KODAMA.chem.sim(Metabolites$SMILE)  
clu=tree.cutting(res,max_nc = 30)
plot(clu$min_nc:clu$max_nc,clu$res.S)
ww=WMCSA(Metabolites$concentration,clu,20)


}
}

