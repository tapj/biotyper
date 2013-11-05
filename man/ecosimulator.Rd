\name{ecosimulator}
\alias{ecosimulator}
\alias{print.ecosimulator}

\title{
simulation phylogenetic or functional composition of hypothetical samples based on models
}
\description{
This function permit to simulate, starting a ecological matrix, phylogenetic or functional composition of hypothetical samples based on uniform, normal and log.normal models.
}
\usage{
ecosimulator(x, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame with samples/sites as columns and species/functions as rows.
}
  \item{r}{
number of simulation
}
}
\details{
For a given set of N feature vectors, each representing the phylogenetic (genus, phylum, etc) or the functional (gene, orthologous group, functional modules, etc) composition of a metagenomic sample, we simulated N hypothetical metagenomic samples containing the same number of features by sampling from a continuum as follows:
template-uniform: each generated feature corresponds to a feature in the real dataset and its values across multiple samples is uniformly distributed between the minimum and maximum abundances of that feature in the real dataset
template-Gaussian (normal): each generated feature corresponds to a feature in the real dataset and its values across multiple samples followed a Gaussian distribution with the same mean and standard deviation as observed in the corresponding real data
template-log.normal: each generated feature corresponds to a feature in the real dataset and its values across multiple samples followed a log.normal distribution with the same mean and standard deviation as observed in the corresponding real data

Generated values are then normalized so that abundances within a sample sum to 1.

}
\value{
a list of class \code{"ecosimulation"}
\item{uniform}{simulated \code{\link{data.frame}} using uniform distribution}
\item{normal}{simulated \code{\link{data.frame}} using normal distribution}
\item{log.normal}{simulated \code{\link{data.frame}} using log.normal distribution}
.
}
\references{
Arumugam et al. Enterotypes of of the human gut microbiome. Nature 2011. doi:10.1038/nature09944
}
\author{
Julien Tap
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{runif}}
\code{\link{rnorm}}
\code{\link{rlnorm}}
}
\examples{

data(Titanium16S)
simul=ecosimulator(Titanium16S, 50) 
simul
class(simul)

}


\keyword{simulation}
