geosub <- function ( gffs , ndffile , setname , samplenames = basename(removeext(gffs)) , cores=length(gffs) , qnorm = TRUE ){
	
	library(parallel)
	
	numgffs<- length(gffs)

	# read in files
	cat("reading in gffs\n")
	gffl<-mclapply(gffs,read.tsv,mc.cores=cores)
	
	cat("reading in ndf\n")
	ndf<-read.tsv(ndffile,header=T)
	ndf$PROBE_ID <- gsub("_RS","RS",ndf$PROBE_ID)


	# sort gffs
	cat("sorting data\n")
	gffl<-mclapply(gffl, function(x) x[order(x$V9),] , mc.cores = cores )
	
	# test if all identical
	cat("testing if all gffs have identical probes\n")
	ids<-gffl[[1]]$V9
	tests <- unlist(mclapply(gffl, function(x) identical(x$V9,ids) , mc.cores = cores ))
	if(all(tests)==FALSE){ stop ("gffs have different probes")}

	probenames<-gsub("_RS","RS",remove.prefix(remove.suffix(ids,";count"),"probe_id="))
	print(head(probenames))
	print(tail(probenames))
	print(length(probenames))
	print(nrow(ndf))
	# combine gffs in single table
	scores<-as.data.frame(mclapply(gffl,"[",6,mc.cores=cores))
	
	mat <- as.data.frame( cbind ( gffl[[1]][,1] , gffl[[1]][,4] , gffl[[1]][,5] , probenames , scores ) , stringsAsFactors=FALSE )

	rm(gffl)

	colnames(mat) <- c("chr","start","stop","PROBE_ID" , samplenames)

	nummatrows <- nrow(mat)
	print(nummatrows)

	cat("merging ndf with gffs\n")
	mat <- merge(mat,ndf,by="PROBE_ID",all.x=TRUE,all.y=FALSE)
	print(nrow(mat))
	mat$start <- as.numeric(mat$start)
	mat$stop <- as.numeric(mat$stop)

	#if(nrow(mat) != nummatrows) { stop("lost rows during merge") }

	mat <- as.data.frame( cbind( mat$PROBE_ID, mat$PROBE_SEQUENCE, mat$chr, mat$start, mat$stop-mat$start, mat$stop, mat[,4+(1:numgffs)]), stringsAsFactors=FALSE)

	colnames(mat) <- c( "ID_REF", "SEQUENCE", "Chromosome", "RANGE_START", "LENGTH", "RANGE_END", samplenames )

	# quantile normalize
	cat("normalizing data\n")
	if(qnorm){
		scores<-mat[,6+(1:numgffs)]
		for(i in 1:numgffs){
			scores[,i] <- scores[order(scores[,i]),i]
		}
		avgdist<-rowMeans(scores)
		avgdist<-avgdist[order(avgdist)]
		for(i in 6+(1:numgffs) ) {
			mat[order(mat[,i]),i]<-avgdist
		}
	}

	# save matrix
	cat("saving matrix\n")
	write.table(mat,file=paste0(setname,"_MATRIX.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)



}