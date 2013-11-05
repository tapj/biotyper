\name{cluster.biotypes}
\alias{cluster.biotypes}
\alias{cluster.biotypes.default}
\alias{cluster.biotypes.dist}

\title{Using Partitioning around medoids to clusters samples for a given k}

\description{
This function permit from a data frame (species as rows and samples as columns
}
\usage{
cluster.biotypes(x, ...)
\method{cluster.biotypes}{default}(x, k)
\method{cluster.biotypes}{default}(x, k)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame or a distance matrix
}
 \item{k}{
an interger number of clusters
}
  \item{\dots}{
not used
}
}
\details{
This function use \code{\link{pam}}.
}
\value{
a vector containing intenger number
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Julien Tap
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pam}}
\code{\link{CH.index}}
}
\examples{
data(Titanium16S)
cluster.biotypes(Titanium16S, k=3)
}

