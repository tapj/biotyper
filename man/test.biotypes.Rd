\name{test.biotypes}
\alias{test.biotypes}
\alias{print.biotypes.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robustness of enterotypes clusters
}
\description{
To assess whether the samples are continuously distributed in the phylogenetic composition space, or whether they predominantly congregate around a few cluster centers
}
\usage{
test.biotypes(biotyper.obj, biotyper.obj.simul)

\method{print}{biotypes.test}(biotypes.test)
\method{plot}{biotypes.test}(biotypes.test)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{biotyper.obj}{
an object of class \code{"biotyper.obj"}, i.e. a biotyper object
}
\item{biotyper.obj.simul}{
an object of class \code{"biotyper.obj.simul"}, i.e. a biotyper object from simulated data
}

\item{biotypes.test}{ a biotypes.test object from \code{test.biotypes}}

}

\value{

A list of element :
\item{silhouette.test}{p.values according uniform, normal, log.normal distribution}
\item{silhouette.simul}{Silhouette simulated coef. according uniform, normal, log.normal distribution simulation}
\item{observed.sil}{Observed silhouette coef. from  \code{biotyper.obj}}


}

\author{
Julien Tap
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

\dontrun{
## Not run:
data(Titanium16S)
Titanium.jsd=dist.JSD(Titanium16S)
ch=CH.index(Titanium16S.jsd)
k=which(ch==max(ch))

ecosimulation=ecosimulator(Titanium16S, 100)
biotyper.obj.simul=biotyper(ecosimulation, k)


biotyper.obj=biotyper(Titanium16S, k, distance.jsd=Titanium16S.jsd)

test=test.biotypes(biotyper.obj, biotyper.obj.simul)
plot(test)

## End (Not run)
}

}

\keyword{test}

