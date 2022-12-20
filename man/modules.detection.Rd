\name{modules.detection}
\alias{modules.detection}
\title{
Detection of modules.
}
\description{
This function calculates the structural similarity between different metabolites and perform hierarchical clustering using the KODAMA algorithm and detect the optimal number of clusters. The procedure is repeated to ensure the robustness of the detection.
}
\usage{


modules.detection  (smiles,
                     repetition=10,
                     k=50,
                     seed=1234,
                     max_nc = 30,
                     dissimilarity.parameters=list(),
                     kodama.matrix.parameters=list(),
                     kodama.visualization.parameters=list(),
                     hclust.parameters=list(method="ward.D"),
                     verbose = TRUE)

}

\arguments{
\item{smiles}{
A list of smile notations for the study metabolites dataset.
}
\item{repetition}{
The number of time the KODAMA analysis is repeated. 
}
\item{k}{
A number of components of multidimensional scaling.
}
\item{seed}{
Seed for the generation of random numbers.
}
\item{max_nc}{
Maximum number of clusters.
}
\item{dissimilarity.parameters}{
Optional parameters for  \code{\link{chemical.dissimilarity}} function.
}
\item{kodama.matrix.parameters}{
Optional parameters for  \code{\link{KODAMA.matrix}} function.
}
\item{kodama.visualization.parameters}{
Optional parameters for  \code{\link{KODAMA.visualization}} function.
}
\item{hclust.parameters}{
Optional parameters for  \code{\link{hclust}} function.
}
\item{verbose}{
If verbose is TRUE, it displays the progress for each iteration.
}
}

\value{
  A list contains all results of KODAMA chemical similarity analysis and hierarchical clustering.

}


\seealso{
  \code{\link{KODAMA.matrix}}
}

\examples{
\donttest{
data(Metabolites)

res=modules.detection(Metabolites$SMILE) 

}
}

