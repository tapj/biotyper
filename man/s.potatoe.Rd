\name{s.potatoe}
\alias{s.potatoe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of factorial maps with represnetation of point classes
}
\description{
performs the scatter diagrams with representation of point classes, using xpline to shape graphical representation.
}
\usage{
s.potatoe(dfxy, fac, xax = 1, yax = 2, col.border = seq(1:length(levels(fac))) + 1, col.fill = seq(1:length(levels(fac))) + 1, shape = -0.5, open = "FALSE", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfxy}{
a data frame containing the two columns for the axes
}
  \item{fac}{
a factor partitioning the rows of the data frame in classes
}
  \item{xax}{
the column number of x in dfxy
}
  \item{yax}{
the column number of y in dfxy
}
  \item{col.border}{
a vector of colors used to draw each class in a different border color
}
  \item{col.fill}{
a vector of colors used to draw each class in a different filled color
}
  \item{shape}{
a number to define the shape. see \code{\link{xpline}}
}
  \item{open}{
see \code{\link{xpline}}
}
  \item{\dots}{
not used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The matched call.
}
\references{
http://pbil.univ-lyon1.fr/R/enseignement.html
}
\author{
Daniel Chessel
Julien Tap
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{s.class}}
}
\examples{
data(Titanium16S)
Titanium16S.biotypes=biotyper(Titanium16S, k=3, manalysis=TRUE)
plot(Titanium16S.biotypes$dPCoA$l1[1:2])
s.potatoe(Titanium16S.biotypes$dPCoA$l1, fac=factor(Titanium16S.biotypes$biotypes))
}

