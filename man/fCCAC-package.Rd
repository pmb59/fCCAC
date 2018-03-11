\name{fCCAC-package}
\alias{fCCAC-package}
\alias{fCCAC}
\docType{package}

\title{ functional Canonical Correlation Analysis to evaluate Covariance between nucleic acid sequencing datasets }
\description{ An application of functional canonical correlation analysis  to assess covariance of nucleic acid sequencing datasets such as  chromatin immunoprecipitation followed by deep sequencing (ChIP-seq).}

\details{
\tabular{ll}{
Package: \tab fCCAC\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2016-09-09\cr
License: \tab Artistic-2.0\cr 
LazyLoad: \tab yes\cr
}
}

\author{
Pedro Madrigal,

Maintainer: Pedro Madrigal \email{dnaseiseq@gmail.com}

}
\references{ Madrigal P (2016) fCCAC: functional canonical correlation analysis to evaluate covariance between nucleic acid sequencing datasets. Bioinformatics: http://doi.org/10.1093/bioinformatics/btw724.  }

\examples{

## hg19. chr21:40000000-48129895 H3K4me3 data from Bertero et al. (2015)
if (.Platform$OS.type == "unix") {
  
owd <- setwd(tempdir())

bigwig1 <- "chr21_H3K4me3_1.bw"
bigwig2 <- "chr21_H3K4me3_2.bw"
bigwig3 <- "chr21_H3K4me3_3.bw"
peakFile <- "chr21_merged_ACT_K4.bed"
labels <- c( "H3K4me3", "H3K4me3","H3K4me3" )
ti <- "H3K4me3 peaks"

r1 <- system.file("extdata", bigwig1,  package="fCCAC",mustWork = TRUE)
r2 <- system.file("extdata", bigwig2,  package="fCCAC",mustWork = TRUE)
r3 <- system.file("extdata", bigwig3,  package="fCCAC",mustWork = TRUE)
r4 <- system.file("extdata", peakFile, package="fCCAC",mustWork = TRUE)

fc <- fccac(bar=NULL, main=ti, peaks=r4, bigwigs=c(r1,r2,r3), labels=labels, splines=15, nbins=100, ncan=15) 

head(fc)



setwd(owd)
}
}
