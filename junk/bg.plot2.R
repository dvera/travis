bgPlot <-
function( bgfiles, regions , allsameylims=TRUE , chromosomes=FALSE , lspan = 0 , overlay=TRUE, linewidth=2, ylabel="score" , ylims=c(NA,NA) , plotcolors=rainbow(length(bgfiles)) , legendnames = basename(removeext(bgfiles)) , threads=getOption("threads",1L) ){

	numbgs <- length(bgfiles)

	bgl <- bgRead(bgfiles,threads=threads,scores=FALSE)

	# if(chromosomes){
	# 	chroms<-mixedsort(unique(unlist(lapply(bgl,"[",1))))
	# 	numchroms<-length(chroms)
	# }

	numregions <- length(regions)

	#all <- do.call(rbind,bgl)
	#maxlength <- max(unlist(lapply(bgl,"[",3)))

	scores <- unlist(lapply(bgl,"[",4))
	maxscore <- max(scores)
	minscore <- min(scores)
	rm(scores)


	#par(mfrow=c((numchroms+1),1),mar=c(0,3,0,2),oma=c(3,5,3,3))

	plot(0,type='n',axes=F,ann=F)
	legend("topright",legend=legendnames,lwd=3,col=plotcolors,horiz=T)
	for(i in 1:numregions){
		curchroms<-lapply(1:numbgs, function(z) bgl[[z]][which(bgl[[z]][,1]==chroms[i]),])
		chromsize<-max(curchroms[[1]][,3])
		allscores<-unlist(lapply(1:numbgs, function(z) curchroms[[z]][,4] ) )
		autoylims<-c(min(allscores,na.rm=TRUE),max(allscores,na.rm=TRUE))
		rm(allscores)

		if(allsameylims & is.na(ylims[1])){
			ylims[1] <- minscore
		}
		if(allsameylims==FALSE & is.na(ylims[1])){
			ylims[1] <- autoylims[1]
		}
		if(allsameylims & is.na(ylims[2])){
			ylims[2] <- maxscore
		}
		if(allsameylims==FALSE & is.na(ylims[2])){
			ylims[2] <- autoylims[2]
		}


		if(overlay==TRUE){
			plot(
				0,
				type="n",
				ylim=ylims,
				xlim=c(0,maxlength),
				ylab=ylabel,
				xlab=if(i!=numchroms){xaxt='n'},
			)


			#if(i!=numchroms){axis(side=1,labels=F)}

			mtext(chroms[i],side=4,las=2,line=1)

			for(j in 1:numbgs){
				if(lspan != 0){
					#plotscores <- loess(curchroms[[j]][,4]~curchroms[[j]][,2],span=lspan*(chromsize/maxlength))$fitted
					plotscores <- loess(curchroms[[j]][,4]~curchroms[[j]][,2],span=lspan*(maxlength/chromsize))$fitted
				} else{
					plotscores <- curchroms[[j]][,4]
				}
				lines(curchroms[[j]][,2],plotscores,lwd=linewidth, col=plotcolors[j])
			}

			grid(nx=NULL,ny=NA,col="grey70")
			abline(h=0,v=seq(0,200,10),col="grey70")

			#if(i==numchroms){paste0("Chromosome coordinate (1 x 10^",ylimordermag," bp)")}

		}

		mtext(ylabel,side=2,outer=T,line=2)


		if(overlay==FALSE){
			par(mfrow=c(numbgs,1),oma=c(4,2,4,0),mar=c(0,4,0.75,10))
			for(j in 1:numbgs){
				plot(
					curchroms[[j]][,2],
					curchroms[[j]][,4],
					type="h",
					ylim=if(ylims[1]=="auto"){autoylims} else{ylims},
					xlim=c(0,chromsize),
					ylab="",
					xlab=if(j==numbgs){"chromosome coordinate (bp)"} else{""},
					xaxt=if(j==numbgs){'s'} else{'n'},
					col=plotcolors[j]

				)
				axis(side=1,labels=F)
				grid(lty="solid",ny=NA)
				abline(h=0,lwd=1,col="black")
				mtext(legendnames[j],side=4,las=1)
			}
			mtext(ylabel,side=2,outer=T,cex=1.5)
			mtext(chroms[i],side=3,outer=T,cex=1.5)
		}
	}
	dev.off()
}
