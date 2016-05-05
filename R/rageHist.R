rageHist <-
function( scorelist , vlines=0, hlines=0, cdf = FALSE , dens=TRUE , drawlegend = TRUE , legendpos="topright" , fraction = FALSE , brks = 50 , reverse = FALSE , xlims=NULL, ylims=NULL, plotcolors=rainbow(length(scorelist)), legendnames=names(scorelist), threads=getOption("threads",1L) , printn=TRUE ) {

	numbeds<-length(scorelist)

	scores <- scorelist

	if(is.null(names(scores))){ names(scores) <- 1:length(scores) }

	if(printn){
		numfrags<-unlist(lapply(scores, length))
	} else{
		numfrags=NULL
	}
	# set x axis lims
	if(!is.null(xlims)){

		if( grepl( "%" , xlims[1] ) ) {
			xlims[1] <- as.numeric(min(unlist(lapply(scores,quantile, probs=as.numeric(gsub("%","",xlims[1] ) ) / 100 , na.rm=T ))))
		}

		if( grepl( "%" , xlims[2] ) ) {
			xlims[2] <- as.numeric(max(unlist(lapply(scores,quantile, probs=as.numeric(gsub("%","",xlims[2] ) ) / 100 , na.rm=T ))))
		}

	} else{

		xlims=c(min(unlist(scores)), max(unlist(scores)))

	}

	xlims <- as.numeric( xlims )

	if(cdf==FALSE){
		if(dens){
			densitylist<-mclapply( scores, density, from=xlims[1], to=xlims[2], mc.cores=threads )
			xl<-lapply( lapply( densitylist, "[" , 1), unlist )
			yl<-lapply( lapply( densitylist, "[" , 2), unlist )
			ymax<-max(unlist(yl))
			ylabel = "density"

		} else{
			scores <- mclapply( 1:numbeds, function(x) scores[[x]][which(scores[[x]] >= xlims[1] & scores[[x]] <= xlims[2] )] )
			densitylist<-mclapply( scores, hist, plot=FALSE, breaks=brks )
			xl<-lapply( lapply( densitylist, "[" , 4), unlist )
			yl<-lapply( lapply( densitylist, "[" , 2), unlist )
			ylabel = "count"
			if(fraction){
				yl <- lapply( 1:numbeds, function(x) yl[[x]]/length(scores[[x]] ) )
				ylabel="proportion"
			}
			ymax<-max(unlist(yl))

		}

	} else{
		xl <- lapply( scores, sort )
		yl <- lapply( 1:numbeds, function(x) length(xl[[x]]) - ( 1:length(xl[[x]]) ) )
		ylabel = "cumulative sum"
		if(reverse){
			yl <- lapply( 1:numbeds, function(x) length(xl[[x]]) - yl[[x]] )
		}
		if(fraction){
			yl <- lapply( 1:numbeds, function(x) yl[[x]]/length(xl[[x]]) )
			ylabel = "cumulative proportion"
		}

		ymax<-max(unlist(yl))

	}

	if(is.null(ylims)){
		ylimits <- c(0,ymax)
	} else{
		ylimits <- ylims
	}

	plot(0,type="n",ylim=ylimits, xlim=xlims, xlab="value", ylab=ylabel)

	for(i in 1:numbeds){
		lines(xl[[i]],yl[[i]],col=plotcolors[i],lwd=2)
	}

	abline(v=vlines, col="grey80")
	abline(h=hlines, col="grey80")

	if(drawlegend){
		legend(legendpos,legend=paste(legendnames,", n=",if(printn){numfrags}), col=plotcolors, lwd=3)
	}

}
