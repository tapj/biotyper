\name{biotyper}
\alias{biotyper}
\alias{biotyper.data.frame}
\alias{biotyper.ecosimulation}
\alias{plot.biotyper}
\alias{print.biotyper.obj.simul}
\alias{summary.biotyper.obj.simul}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biotyping your commununity}
\description{This function permit to cluster your sites/samples (columns) in k biotypes according distribution of features (species or functions in rows), calculate cluster and clasification index, provide multivariate analysis as PCA, between class and dPCoA, and plot everything. The same function take in charge ecosimulation object.}

\usage{
biotyper(x, ...)

\method{biotyper}{data.frame}(obs, k=3 , distance.jsd=NULL, cluster=NULL, manalysis=FALSE,  plot=FALSE, nf=3, no.unassigned=TRUE,...)
\method{biotyper}{ecosimulation}(obs, k , ...)
\method{plot}{biotyper}(biotyper.obj, xax=1, yax=2, potatoes=FALSE, ...)
\method{print}{biotyper.obj.simul}(biotyper.obj.simul)
\method{summary}{biotyper.obj.simul}(biotyper.obj.simul)

}

\arguments{

  \item{obs}{ A data.frame with sites as columns and species/function as rows.}
  \item{k}{number of clusters}
  \item{distance.jsd}{a distance matrix of sites.if NULL a jsd distance matrix will be created using dist.JSD  }
  \item{cluster}{a numeric or character vector to force your cluster}
  \item{manalysis}{logical. if FALSE, no multivariate analysis are done}
  \item{plot}{logical. if FALSE, nothing is plotted}
  \item{nf}{Numbers of component, if you do not know keep nf=3 by default.}
  \item{no.unassigned}{logical. Do you have an unassigned features (on the first row) ? If TRUE, multivariate analysis will remove it during PCA analysis.}
  \item{\dots}{not used}
  \item{biotyper.obj}{an object of class \code{"biotyper.obj"}, i.e. a biotyper object}
  \item{biotyper.obj.simul}{an object of class \code{"biotyper.obj.simul"}, i.e. a biotyper object}
  \item{potatoes}{logical. if TRUE, \code{\link{s.potatoes}} will be used to represent biotypes}
  \item{xax}{the component number for the x axis}
  \item{yax}{the component number for the y axis}


 }



\details{
Have a look to Arumugam et al.Nature.2011.
}
\value{

Acoording the class of object (i.e a \code{data.frame} or \code{ecosimulation}), an object of class \code{biotyper.obj} or \code{biotyper.obj.simul}, basically including a list including elements
\item{biotypes}{input or generated biotypes using \code{\link{pamk}} function from code{\link[fpc:fpc-package]{fpc}}}
\item{obs.silhouette}{silhouette coefficient}
\item{PCA}{Principal component analysis using \code{\link[ade4:ade4-package]{ade4}}}
\item{BET}{Between class analysis using \code{\link[ade4:ade4-package]{ade4}}}
\item{dPCoA}{Double principle coordinate analysis using \code{\link[ade4:ade4-package]{ade4}}}


}
\references{
Arumugam et al. Enterotypes of of the human gut microbiome. Nature 2011. doi:10.1038/nature09944
}
\author{
Julien Tap
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
objects to See Also as \code{\link{test.biotypes}} , \code{\link{CH.index}} , \code{\link{ecosimulator}}
}
\examples{
data(Titanium16S)
Titanium16S.jsd=dist.JSD(Titanium16S)
ch=CH.index(Titanium16S,Titanium16S.jsd,k=1:10)
k=which(ch==max(ch))
biotyper(Titanium16S, k, distance.jsd=Titanium16S.jsd, manalysis=TRUE, plot=TRUE)


\dontrun{
## Not run:
data(Titanium16S)
Titanium.jsd=dist.JSD(Titanium16S)
ch=CH.index(Titanium16S,Titanium16S.jsd,k=1:10)
k=which(ch==max(ch))
ecosimulation=ecosimulator(Titanium16S, 100)
biotyper.obj.simul=biotyper(ecosimulation, k)
summary(biotyper.obj.simul)
## End (Not run)
}

}


\keyword{clustering}
\keyword{biotypes}
\keyword{enterotypes}







