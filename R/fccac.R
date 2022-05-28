
fccac <- function(peaks, bigwigs, labels, splines=10, nbins=100, ncan=5 , tf=c(), main="", bar=NULL, outFiles=FALSE ){    
		
	# Read genomic regions in BED format with Bioconductor package genomation
	print("Reading peaks...")
	peaks <- readGeneric(file=peaks, chr=1, start=2, end=3)  #[1:100]  #readBed
	  
	if (ncan > splines | ncan > length(peaks) ){ print("ncan should not be higher than the number of splines or peaks. Please lower the value of ncan.")   }

	if (ncan <= splines | ncan <= length(peaks)){
		
		print("Starting fCCAC...")
		SPLINES <- splines
		L <- nbins    # bp

		# Modify labels to detect replicates (assume replicates are ordered 1,2,...)
		d <- duplicated(labels)
		new_labels <- rep("newLabel", length(labels) ) 
		for (l in seq(from=1, to=length(labels), by=1)  ) {
			if (d[l] == FALSE) {  counter <- 1; new_labels[l]  <- paste(labels[l], 1, sep="_Rep") }
			if (d[l] == TRUE)  {  counter <- counter +1;  new_labels[l]  <- paste(labels[l], counter, sep="_Rep")  } 	
		}
	
		# Create B-spline basis system (cubic, order 4) of the signals before FDA:
		bspl <- create.bspline.basis(rangeval=c(-L/2,L/2), nbasis=SPLINES, norder=4)   # Cubic B-splines
		argvalsBS <- seq(from= (-L/2), to=(L/2), by=1)      #(-L/2):(L/2)
		
		# Read bigWig Files 
		fdaData <- lapply(seq(from=1, to=length(bigwigs), by=1), function(x) matrix(NA, nrow=length(peaks) , ncol=L+1 ))

		for (i in seq(from=1, to=length(bigwigs), by=1)   ) {
			print(paste(c("Reading bigWig file...",i,"/",length(bigwigs)),collapse="")  )

			fdamatrix <- matrix(0.0, ncol=L+1, nrow= length(peaks) ) 
			fdamatrix  <- ScoreMatrixBin(target = bigwigs[i], bin.num = L+1, windows = peaks, type="bigWig",rpm=FALSE, strand.aware = TRUE, bin.op="max" )
			fdaData[[i]] <- Data2fd(y=t(fdamatrix), argvals= argvalsBS, basisobj=bspl) 			
		}
		rm(fdamatrix)
		# length(fdaData)

		co <- combn(x=c(new_labels), m=2)  	# all possible pairwise combinations

		# Select Sample of Interest (if any)
		if(length(tf)==1){
    			idx <- unique(col(co)[grepl(tf, co)])
    			co <- co[, idx] 
  		}
		
		# Prepare Data for functional CCA
		scc <- rep(0, ncan*ncol(co) )  		# squared Canon. Corr.
		sccM <- rep(0, ncan*ncol(co) ) 		# MAX squared Canon. Corr.
		pair <- rep("pair", ncan*ncol(co) ) 	# labeling
		Spair <- rep("Spair", ncol(co) ) 	# labeling for S
			
		w <-  1/ ( seq(from=1, to=ncan, by=1)  )
		Ma <-  sum(w) 				# maximum possible value for w
		S <- rep(0, ncol(co) )  		# weighted sums
		x <- lapply(seq(from=1, to=ncol(co), by=1)  , function(x) NA) 	# list to store output of cca.fd

		for (i in seq(from=1, to=ncol(co), by=1)  ) {
			print(paste(c("Performing fCCA in pair...",i,"/",ncol(co)) ,collapse="") )
		
			file1 <- which(new_labels==co[1,i])
			file2 <- which(new_labels==co[2,i])
		
			print(paste(c(new_labels[file1],"...vs...",new_labels[file2]) , collapse="") )
	
			x[[i]] = cca.fd(fdobj1=fdaData[[file1]], fdobj2=fdaData[[file2]], ncan = ncan, ccafdPar1=fdPar(bspl, 0, 0), ccafdPar2=fdPar(bspl,0, 0), centerfns=TRUE)
	
		}
		# NCAN calculated in 'cca.fd' is always the number of splines
		# plot.cca.fd(x)

		for (i in seq(from=1, to=ncol(co), by=1)  ) {
	
			posit <- (1+(ncan*(i-1))):(ncan*(i-1+1) )
			
			scc[posit] <-  x[[i]]$ccacorr[ seq(from=1, to=ncan, by=1) ]                       # all canonical correlations
			sccM[posit] <-  rep( max (x[[i]]$ccacorr[ seq(from=1, to=ncan, by=1) ]), ncan)    # max of all canonical correlations

			file1 <- which(new_labels==co[1,i])
			file2 <- which(new_labels==co[2,i])
			pair[posit] <-  rep( paste(new_labels[c(file1,file2)]  , collapse="_vs_") ,ncan)  
			Spair[i] <- paste(new_labels[c(file1,file2)], collapse="_vs_" ) 

			# calculate weigthed sum
			S[i] <- sum(w* x[[i]]$ccacorr[ seq(from=1, to=ncan, by=1) ]    )
			
			rm(posit, file1, file2)

		}
		# if ( is.null(bar) == TRUE ) {bar <- ncol(co) }  ##barplot(100* S/Ma, ylim= c(0,100)  )

		# Colormap
		colfunc <- colorRampPalette(brewer.pal(10, "Spectral"))

		# Plots
		# ggplot2
		ggData <- data.frame(scc = scc, pair = pair, variables = rep( seq(from=1, to=ncan, by=1) , ncol(co)), sccM = sccM)	
		ggData <- transform(ggData,  pair = reorder(pair, sccM))
		p1 <- ggplot(ggData, aes(x=variables, y=scc, group=pair))  + geom_line(aes(colour = pair)) + ylim(0,1)  + theme_bw() + theme(panel.border = element_rect( colour = "black")) +   theme(legend.position='none', plot.title = element_text(lineheight=.8, face="bold")) + xlab("Canonical variable (k)") + ylab( expression(Squared~canonical~correlation~(R[k] ^{2}) ) )  +   geom_point(size=0.9, shape=5, aes(colour=pair)) + scale_x_continuous(breaks=seq(from=1, to=ncan, by=1)  ) + scale_colour_manual(values = colfunc(ncol(co))) + ggtitle(main)  
		# select 1st canonical correlation
		ggData <- subset(ggData, variables==1 )
		# sort and assign colors
		ggData <- ggData[sort(ggData$sccM, index.return=TRUE, decreasing=TRUE)[[2]],]
		ggData$CL <- rev( colfunc(ncol(co))  )

		ggData2 <- data.frame(S = 100* (S/Ma), pair = Spair  )   
		#ggData2 <- transform(ggData2,  pair = reorder(pair, S))
		ggData2 <- ggData2[sort(ggData2$S, index.return=TRUE, decreasing=TRUE)[[2]],]
		ggData3 <- merge(x=ggData2, y=ggData, by='pair',  sort = FALSE)
		#ggData3 <- ggData3[sort(ggData3$S, index.return=T, decreasing=T)[[2]],]
		ggDataTXT <- ggData3
			
		INB = is.null(bar)
		if ( INB == TRUE ) {   bar <- ncol(co); CHOSEN <- seq(from=1, to=bar, by=1)  ;  ggData3 <- ggData3[CHOSEN,]  }
		if ( INB != TRUE ) {  CHOSEN<-c(seq(from=1, to=bar[1], by=1),  (length(ggData3$CL)-bar[2]+1):length(ggData3$CL)); ggData3 <- ggData3[  CHOSEN , ]	}
		
		p2 <-  ggplot(ggData3, aes(x = reorder(factor(pair),S), y = S)) + geom_bar(stat = "identity",width=.6, fill=reorder(rev(as.character(ggData3$CL)) ,ggData3$S) )+ coord_flip() + ylim(0,100)  + theme_bw() + theme(panel.border = element_rect( colour = "black")) +  theme(legend.position='none', text = element_text(size=10)   )  + ylab("F(%)") + xlab("") + geom_hline(yintercept = 100, colour="red", linetype = "longdash")
		
		colnames(ggDataTXT) <- c("samples","F","squ_can_corr_k_1","k","squ_can_corr_k_1","color"  )
		fccac_out <- ggDataTXT[,c(1,2,3,4,6)]
		
		if (outFiles == TRUE){
			print("Saving fCCAC.pdf...")
			pdf("fCCAC.pdf", height=6, width=3.5)	
			multiplot(p1, p2, cols=1)
			dev.off()
			print("Saving fCCAC.txt...")
			write.table(x=fccac_out , file = "fCCAC.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
			print("Done...")
		}

		if (outFiles == FALSE){
			multiplot(p1, p2, cols=1)
		}		
		return( fccac_out )
	}	
}
