\name{CH.index}
\alias{CH.index}

\title{Evaluate the optimal number of cluster using CH index}
\description{This function permit from a ecological data frame and a distance matrix between samples to estimated CH index for a range of k clusters}
\usage{
CH.index(dataframe, distance.JSD, kvector = 1:20, clusterSim=FALSE)
}

\arguments{
  \item{dataframe}{
A data frame with species/functions as rows and samples/sites as columns
}
  \item{distance.JSD}{
a distance matrix between samples/sites
}
  \item{kvector}{
a numerical vector of k cluster, by default 1:20
}
  \item{clusterSim}{
logical: if TRUE clusterSim::indexG1 will be used instead of pamk
}
}
\details{
This functions allow you to estimate the optimal number of biotypes in your dataset taking the k corresponding to the highest CH index.
}
\value{
\item{vector}{a vector a CH index for each k}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Julien Tap
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{index.G1}}
\code{\link{cluster.biotypes}}

}
\examples{
data(Titanium16)
CH.index(Titanium16S, k=1:10) -> ch
plot(ch, xlab=" k clusters", ylab="CH index")
abline(v=which(ch==max(ch)), lwd=2, col="red") #optimal k

}

\keyword{CH}

