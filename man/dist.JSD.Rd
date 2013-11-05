\name{dist.JSD}
\alias{dist.JSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jensen Shannon Distance Matrix Computation
}
\description{
This function computes and returns the distance matrix computed by using JSD distance measure to compute the distances between the column of a data matrix. 
}
\usage{
dist.JSD(inMatrix, pseudocount = 10^(round(log10(min(as.matrix(inMatrix)[as.matrix(inMatrix) > 0])), 0) - 1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inMatrix}{
a data matrix
}
  \item{pseudocount}{
a pseudocount used to replace zero count
}
  \item{\dots}{
not used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{\link{dist.JSD}} returns an object of class \code{"dist"}. 
}

\author{
Julien Tap
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

\code{\link{dist}}

}
\examples{
data(Titanium16S)
Titanium16S.jsd=dist.JSD(Titanium16S)
class(Titanium16S.jsd)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distance}
\keyword{JSD}
