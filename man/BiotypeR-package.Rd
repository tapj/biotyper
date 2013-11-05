\name{BiotypeR-package}
\alias{BiotypeR-package}
\alias{BiotypeR}
\docType{package}
\title{ BiotypeR : an R package to biotype a community }
\description{
This package provide numerous functions for biotyping your community dataset based on clustering and classification techniques.
}

\details{
\tabular{ll}{
Package: \tab BiotypeR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-07-12\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}

}
\author{
Julien Tap
Maintainer: Julien Tap <tap@embl.de>
}
\references{
Arumugam et al. Enterotypes of of the human gut microbiome. Nature 2011. doi:10.1038/nature09944
}
%% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%% ~~ the R documentation directory ~~
\keyword{ecology}
\keyword{metagenomics}
\keyword{enterotypes}

\seealso{
Optional links to other man pages, e.g.
\code{\link[<ade4>:<ade4>-package]{<ade4>}} ~~
\code{\link[<fpc>:<fpc>-package]{<fpc>}} ~~
}
\examples{
data(Titanium16S)
Titanium16S.biotypes=biotyper(Titanium16S, k=3, manalysis=TRUE)
plot(Titanium16S.biotypes)
}
