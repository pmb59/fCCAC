\name{fccac}
\alias{fccac}

\title{ functional Canonical Correlation Analysis to evaluate Covariance between nucleic acid sequencing datasets }
\description{ functional Canonical Correlation Analysis to evaluate Covariance between nucleic acid sequencing datasets. }

\usage{
fccac(peaks, bigwigs, labels, splines=10, nbins=100, ncan=5 , tf=c(), main="", 
bar=NULL, outFiles=FALSE )
}

\arguments{
  \item{peaks}{ BED file. Column 1: chr, Column 2: start, Column 3: end (REQUIRED). }
  \item{bigwigs}{ A vector of characters containing the path to bigwigs files. Replicates of the same samples should be entered consecutive one another (REQUIRED).  }
  \item{labels}{ IDs for each sample. Replicates should have the same label and be ordered (vector of characters, REQUIRED) }
  \item{splines}{ Number of cubic B-splines used to smooth the data and to estimate the canonical variate weight functions (default: 15) }
  \item{nbins}{ Integer value representing the number of bins that should be used for each window (default: 100) }
  \item{ncan}{ Number of canonical components to report in the results. It cannot be higher than number of splines or the number of peaks (default: 15) }
  \item{tf}{ Plot results involving only this TF or TF-replicate (character). Eg., "SOX2" or "SOX2{\_}Rep1" (default: empty vector. plot all) }
  \item{main}{ Title of the plot generated (default: no title) }
  \item{bar}{ In the barplot, plot only first bar[1] and last bar[2] interactions after ranking by F-value (default: NULL, plots all the combinations).  }
  \item{outFiles}{If TRUE, the function writes two files in the working directory, fCCAC.pdf and fCCAC.txt (tabulated text-file with results). (default: FALSE) }
  
}

\details{ Detailed information about the methodology can be found in Madrigal (2016).}

\value{ The function reports a dataframe with the following columns: pairwise samples, F value, k (order of the first canonical correlation), and value of the first canonical correlation. }
  



\references{ Madrigal P (2016) fCCAC: functional canonical correlation analysis to evaluate covariance between nucleic acid sequencing datasets. Bioinformatics: http://doi.org/10.1093/bioinformatics/btw724.  }

\author{ Pedro Madrigal, \email{pm12@sanger.ac.uk} }

\seealso{ \code{\link{fCCAC-package}} }

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

setwd(owd)
}
}
