bed.spacing <-
function( beds, sizerange=c(1,500), distrange=c(1,1000), numfrags="all", legendnames=basename(removeext(beds)),plotcolors=rainbow(length(beds)),reference="center", cores="max", readsperblock=1000 ){
	library(parallel)
	if(reference %in% c("center","end") == FALSE){stop("not a valid reference point, use 'center' or 'end'")}
	if(cores=="max"){cores=detectCores()-1}
	numbeds<-length(beds)
	bednames<-basename(removeext(beds))
	
	
	beds<-unlist(mclapply(1:numbeds, function(i){
		
		if(is.infinite(sizerange[2]) == FALSE & sizerange[1]>1){
			beds[i]<-bed.parselengths(beds[i],brks=sizerange)
		}

		if(numfrags != "all"){
			cat(bednames[i],": grabbing",numfrags,"random fragments\n")
			newname<-paste(bednames[i],"_",numfrags,"random.bed",sep="")
			system(paste("randomLines",beds[i],numfrags,newname))
			beds[i]<-newname
		}

		if(reference=="center"){
			beds[i]<-bed.centers(beds[i])
		}
		beds[i]
	},mc.cores=detectCores()))
	
	
	bedcoords<-mclapply(1:numbeds, function(i){
		cat(bednames[i],"loading read coordinates\n")
		as.numeric(readLines(pipe(paste("cut -f 2",beds[i]))))
	},mc.cores=detectCores())
	
	dbedcoords<-mclapply(1:numbeds, function(i){
		cat(bednames[i],"loading duplicate-start read coordinates\n")
		as.numeric(readLines(pipe(paste("cut -f 1,2",beds[i],"| uniq -D | cut -f 2"))))
	},mc.cores=detectCores())
	
	distances<-lapply(1:numbeds, function(i){
		coords<-bedcoords[[i]]
		cat(bednames[i],"measuring distances for all reads\n")
		numblocks<-floor(length(coords)/1000)
		dh<-as.data.frame(mclapply(1:numblocks,function(x){
			d<-dist(coords[((((x-1)*1000):(x*1000))+1)])
			d<-na.omit(as.vector(d[d<=distrange[2] & d>=distrange[1] ]))
			h<-hist(d,breaks=distrange[1]:(distrange[2]+1),plot=FALSE)$counts/1000
			return(h)
		},mc.cores=cores))
		rowMeans(dh)
		
	})
	dmax<-quantile(unlist(distances),probs=0.98)
	
	ddistances<-lapply(1:numbeds, function(i){
		coords<-dbedcoords[[i]]
		cat(bednames[i],"measuring distances for all duplicate-start reads\n")
		numblocks<-floor(length(coords)/readsperblock)
		dh<-as.data.frame(mclapply(1:numblocks,function(x){
			d<-dist(coords[((((x-1)*readsperblock):(x*readsperblock))+1)])
			d<-na.omit(as.vector(d[d<=distrange[2] & d>=distrange[1] ]))
			h<-hist(d,breaks=distrange[1]:(distrange[2]+1),plot=FALSE)$counts/readsperblock
			return(h)
		},mc.cores=cores))
		rowMeans(dh)
		
	})
	ddmax<-quantile(unlist(ddistances),probs=0.98)

	cat("plotting all-read distogram\n")
	plot(0,type="n",xlim=c(distrange[1],distrange[2]),ylim=c(0,dmax),xlab="Distance (bp)",ylab="count", main="phasogram for all reads")
	abline(v=c((0:15)*10),lwd=1,col="grey80")
	for(i in 1:numbeds){
		lines(distrange[1]:distrange[2],distances[[i]],col=plotcolors[i],lwd=2)
	}
	legend("topright",legend=legendnames, col=plotcolors, lwd=3)
	
	X11()
	
	cat("plotting duplicate-start read distogram\n")
	plot(0,type="n",xlim=c(distrange[1],distrange[2]),ylim=c(0,ddmax),xlab="Distance (bp)",ylab="count", main="phasogram for duplicate-start reads")
	abline(v=c((0:15)*10),lwd=1,col="grey80")
	for(i in 1:numbeds){
		lines(distrange[1]:distrange[2],ddistances[[i]],col=plotcolors[i],lwd=2)
	}
	legend("topright",legend=legendnames, col=plotcolors, lwd=3)
}
