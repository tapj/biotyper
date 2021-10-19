## R package to biotype a community
Package: biotyper

Type: Package

Title: biotyper: an R package to biotype a community

Version: 0.1.3

Depends: R (>= 1.8.0), ade4, fpc, clusterSim

Date: 2020-08-27

Author: Julien Tap, Daniel R. Mende

Maintainer: Julien Tap <jtap@jouy.inra.fr>, Daniel R. Mende <mende@embl.de>

Contributor: John Bouranis @bouranij

Description: This package provide numerous functions for biotyping your community dataset based on clustering and classification techniques.

License: GPL 2

LazyLoad: yes

Short story : package developped for the enterotype paper (Arumugam et al., Nature, 2010) at EMBL.
I am still maintain it at Danone Nutricia Research hoping that it could be useful for other users too. more information with the enterotype tutorial : http://enterotype.embl.de/. A more recent method to cluster enterotypes using reference-based assignments can be found here http://enterotypes.org/

Install : you can install the BiotypeR package from GitHub using `devtools` package.

    require(devtools)
    install_github("tapj/BiotypeR")
    library(BiotypeR)
    data(Titanium16S)
    Titanium16S.biotypes=biotyper.data.frame(Titanium16S, k=3, manalysis=TRUE)
    
  
  
  
  
  

