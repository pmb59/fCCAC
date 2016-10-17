test_fCCAC <- function() {
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

	fc <- fccac(bar=NULL, main="H3K4me3 peaks", peaks=r4, bigwigs=c(r1,r2,r3), labels=labels, splines=15, nbins=100, ncan=15, outFiles=FALSE) 


        #check that if minimal input is not introduced the package recognizes an error situation
	#peakFile
	checkException( fccac(bar=NULL, main="H3K4me3 peaks", peaks=c(), bigwigs=c(r1,r2,r3), labels=labels, splines=15, nbins=100, ncan=15)  )
	#bigwigs
	checkException( fccac(bar=NULL, main="H3K4me3 peaks", peaks=r4, bigwigs=c(), labels=labels, splines=15, nbins=100, ncan=15) )
	#labels
	checkException( fccac(bar=NULL, main="H3K4me3 peaks", peaks=r4, bigwigs=c(r1,r2,r3), labels=c(), splines=15, nbins=100, ncan=15) )
	
	setwd(owd)
		
	}	

}
