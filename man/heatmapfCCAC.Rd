\name{heatmapfCCAC}
\alias{heatmapfCCAC}

\title{ Heatmap of F values obtained by Canonical Correlation Analysis }
\description{ Heatmap of F values obtained by Canonical Correlation Analysis. This function can only be used if all pairwise comparisons were computed previously with the function 'fccac', i.e., using "tf=c()". }

\usage{
heatmapfCCAC(fc)
}

\arguments{
  \item{fc}{ Output of the function 'fccac'. }

}


\value{ Plots a Heatmap of F values using the package 'ComplexHeatmap'. }
  



\references{ Madrigal P (2016) fCCAC: functional canonical correlation analysis to evaluate covariance between nucleic acid sequencing datasets. Bioinformatics: http://doi.org/10.1093/bioinformatics/btw724.  }

\author{ Pedro Madrigal, \email{dnaseiseq@gmail.com} }

\seealso{ \code{\link{fccac}} }

\examples{

## hg19. chr21:40000000-48129895 H3K4me3 data from Bertero et al. (2015)
if (.Platform$OS.type == "unix") {
  
owd <- setwd(tempdir())

bigwig1 <- "chr21_H3K4me3_1.bw"
bigwig2 <- "chr21_H3K4me3_2.bw"
bigwig3 <- "chr21_H3K4me3_3.bw"
peakFile <- "chr21_merged_ACT_K4.bed"
labels <- c( "H3K4me3", "H3K4me3","H3K4me3" )

r1 <- system.file("extdata", bigwig1,  package="fCCAC",mustWork = TRUE)
r2 <- system.file("extdata", bigwig2,  package="fCCAC",mustWork = TRUE)
r3 <- system.file("extdata", bigwig3,  package="fCCAC",mustWork = TRUE)
r4 <- system.file("extdata", peakFile, package="fCCAC",mustWork = TRUE)
ti <- "H3K4me3 peaks"

fc <- fccac(bar=NULL, main=ti, peaks=r4, bigwigs=c(r1,r2,r3), labels=labels, splines=15, nbins=100, ncan=15) 

head(fc)

heatmapfCCAC(fc)

setwd(owd)

}
}
