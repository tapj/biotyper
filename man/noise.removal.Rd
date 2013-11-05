\name{noise.removal}
\alias{noise.removal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove noise species/functions
}
\description{
This functions permit to remove low abundant species or function according a percentage or if they don't belong to the defined top ranking species.
}
\usage{
noise.removal(dataframe, percent = NULL, top = NULL, bysample = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
a ecological data.frame
}
  \item{percent}{
a number between 0-100 %
}
  \item{top}{
a number between 1 to the maximum of species/functions. this functions use this criteria in priority.
}
  \item{bysample}{
Logical. if \code{FALSE} the function will take the average number across samples. 
}
}

\value{
a \code{data.frame}
}

\author{
Julien Tap
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(Titanium16)
noise.removal(Titanium16S, top=10)
}

\keyword{noise}

