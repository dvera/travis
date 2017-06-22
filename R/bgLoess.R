bgLoess <-
function( bgfiles, lspan, chromsizes, threads=getOption("threads",1L), mbg=FALSE ){

	options(scipen=9999)
	# assumes sorted bg
	if(lspan<1){
		stop("lspan must indicate distance in bp")
	}
	if(missing(chromsizes)){
		chromsizes<-getOption("chromsizes",NULL)
		if(is.null(chromsizes)){warning("No chromsizes file defined. Using chromosome coordinate ranges to estimate chromosome sizes")}
		chrom.sizes <- tsvRead(chromsizes)
	} else{
		chrom.sizes <- NULL
	}
	

	bgnames <- basename(removeext(bgfiles))
	numbgs <- length(bgfiles)
	outnames <- paste(bgnames,"_loess",lspan,".bg",sep="")
	chromthreads <- floor(threads/numbgs)
	if(chromthreads==0){chromthreads=1}

	#dump <- mclapply(seq_len(numbgs), function(x){
	dump <- lapply(seq_len(numbgs), function(x){
		curbg <- as.data.frame( read.table ( bgfiles[x], header=mbg, stringsAsFactors=FALSE ) )
		
		if(mbg){
			columnnames=colnames(curbg)
		}

		chroms    <- unique(curbg[,1])
		numchroms <- length(chroms)
		windowsize <- median(curbg[,3]-curbg[,2])
		all=split(curbg,curbg[,1])
		
		#lscores<-mclapply(1:numchroms,function(i){
		lscores<-lapply(1:numchroms,function(i){
			
			cur <- all[[i]]
			
			curchrom <- cur[1,1]
			
			if(!is.null(chrom.sizes)){
				chromi <- grep(curchrom,chrom.sizes[,1])
				if(length(chromi)==0){
					stop(paste("Cannot find",curchrom,"in chromsizes"))
				} else{
					chromlspan <- lspan/chrom.sizes[chromi,2]
				}
			} else{
				chromlspan <- lspan/sum(max(cur[,3])-min(cur[,2]))
			}
			
			
			cura <- as.data.frame(lapply(4:ncol(cur), function(k){
				cur[,k] <- tryCatch({
						loess(cur[,k]~cur[,2],span=chromlspan)$fitted
					},warning = function(war){
						print(paste("warning for file",bgnames[x],"chromosome",cur[1,1],":",war))
						out <- loess(cur[,k]~cur[,2],span=chromlspan)$fitted
						return(out)
					},error = function(err){
						print(paste("smoothing failed for file",bgnames[x],"chromosome",cur[1,1],":",err))
						return(cur[,k])
					}
				)
			}))
			
			cura <- cbind(cur[,1:3],cura)
			colnames(cura) <- colnames(cur)
			return(cura)
		#},mc.cores=chromthreads, mc.preschedule=FALSE)
		})
		smoothbg<-do.call(rbind,lscores)
		write.table(smoothbg,file=outnames[x],col.names=mbg,quote=FALSE, sep="\t", row.names=F)
		rm(curbg)
		gc()
	#}, mc.cores=threads, mc.preschedule=FALSE)
	})
	return(outnames)
}
