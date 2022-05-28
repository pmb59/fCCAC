heatmapfCCAC <- function(fc){
   
   M <- rep("name", nrow(fc)*2)  #store sample names
   k <- 0
   for (i in seq(from=1, to=length(fc$samples), by=1)    ){
      temp_names <- strsplit( as.character(fc$samples[i]), split="_vs_")[[1]] 
      k <- k +1 ; M[k] <- temp_names[1]
      k <- k +1 ; M[k] <- temp_names[2]
   }
   M <- unique(M)
  
   #plot heatmap
   Fv <- matrix(NA, nrow=length(M), ncol=length(M))
   for (i in seq(from=1, to=length(M), by=1)  ){
      for (j in seq(from=1, to=length(M), by=1) ){
         o1 <- which(as.character(fc$samples) == paste(M[i], M[j], sep="_vs_" ) )
         o2 <- which(as.character(fc$samples) == paste(M[j], M[i], sep="_vs_" ) )
         if (length(o1)==1){ Fv[i,j] <- fc[o1,2] }
         if (length(o2)==1){ Fv[i,j] <- fc[o2,2] }
         if (length(o1)==0 & length(o2)==0) {Fv[j,i] <- 100 }
      }
   }
   colnames(Fv) <- M
   rownames(Fv) <- M

   Heatmap(Fv, name='F (%)')
}
