bgPlot <- function( bgFiles , regions=NA , plotcolors=rainbow(length(bgFiles)), legendnames=basename(removeext(bgFiles)) ,  lines=T , steps=F , dots=F , ylims=NA , threads=getOption("threads",1L) , linetypes=1, linewidths=3, lspan=0, flank=0 , ylabel="score" , xline=0 , maxplots=100, connectWithin=0 ) {

  options(scipen=9999)

  numbgs=length(bgFiles)
  if(numbgs>length(linetypes)){linetypes=rep(linetypes,numbgs)}
  if(numbgs>length(linewidths)){linewidths=rep(linewidths,numbgs)}

  if(all(file.exists(regions))){
    region <- read_tsv(regions,col_names=FALSE)
    if(ncol(region)<3){region <- data.frame(region[,1],1,region[,2],stringsAsFactors=FALSE)}
  } else if(grepl("-",regions) & grepl(":",regions)){
    region <- gsub(",","",regions)
    region <- gsub(":","-",region)
    region <- strsplit(region,"-")
    region <- as.data.frame(t(as.data.frame(region,stringsAsFactors=F)),stringsAsFactors=F)
    region[,2] <- as.numeric(region[,2])
    region[,3] <- as.numeric(region[,3])
    #print(region)
  }

  numplots <- nrow(region)
  if(numplots>maxplots){numplots <- maxplots}

  for(r in seq_len(numplots)){
    #print(region[r,3])
    print(paste("region",r,region[r,1],region[r,2]-flank,region[r,3]+flank))
    # fix this to handle no returned scores
    cmdString <- paste0(
      "printf \"%s\\t%s\\t%s\" \"", region[r,1],"\" \"",region[r,2]-flank,"\" \"",region[r,3]+flank,"\" ",
      "| bedtools intersect -u -a ", bgFiles, " -b stdin "
    )


    scores <- cmdRun( cmdString, tsv=TRUE, threads=threads)

    #if(nrow(scores))

    xlims=c( min(unlist(lapply(scores,"[",2))), max(unlist(lapply(scores,"[",2))))
    if(is.na(ylims)){ylimits=c( min(unlist(lapply(scores,"[",4))), max(unlist(lapply(scores,"[",4))))} else{ylimits=ylims}

    plot(0,0,type='n', xlim=xlims, ylim=ylimits, xlab=paste(region[r,1],"coordinate (bp)") , ylab=ylabel , main=paste0(region[r,1],":",region[r,2],"-",region[r,3]) )
    if(is.numeric(xline)){abline(h=xline)}
    for(i in seq_len(numbgs)){

      xvals <-  scores[[i]][,2]
      numvals <- length(xvals)
      xvals2 <- scores[[i]][,3]
      yvals <-  scores[[i]][,4]
      adjacent <- which( xvals[2:length(xvals)] - xvals2[1:(length(xvals2)-1)] <= connectWithin )
			xvals3 <- (xvals+xvals2)/2

      if(lspan>0){
        yvals <- tryCatch({loess(yvals~xvals,span=lspan)$fitted}, error=function(w){yvals} )
      }

      if(dots){
        points(xvals3,yvals,col=plotcolors[i], pch=20)
      }
			if(steps){
				segments(xvals,yvals,xvals2,yvals,col=plotcolors[i], lty=linetypes[i], lwd=linewidths[i])
			}
      if (lines){
        #segments(xvals,yvals,xvals2,yvals,col=plotcolors[i], lty=linetypes[i], lwd=linewidths[i])
        lines(xvals3,yvals,col=plotcolors[i], lty=linetypes[i], lwd=linewidths[i])
      }

      # if(length(adjacent)>0 & linewidths[i]>0){
        #adjacent <- adjacent[adjacent!=numvals]
        #segments(xvals2,yvals[adjacent],xvals2,yvals[adjacent+1],col=plotcolors[i], lty=linetypes[i], lwd=linewidths[i])
      # }



      if(flank !=0){abline(v=region[r,2:3])}
      #segments(xvals,yvals,xvals2, yvals, col=plotcolors[i], lty=linetypes[i], lwd=linewidths[i])
    }
    legend("topright",legend=legendnames, lwd=linewidths+1, cex=0.6, col=plotcolors, lty=linetypes)


    # if(length(unique(unlist(lapply(scores,nrow))))==1){
    #   sss=as.data.frame(lapply(scores,"[",4))
    #   colnames(sss)<-legendnames
    #   sss[sss>ylims[2]]<-ylims[2]
    #   sss[sss<ylims[1]]<-ylims[1]
    #   brks=seq(ylims[1],ylims[2],(ylims[2]-ylims[1])/100)
    #   numbreaks=length(brks)
    #   heatmap.2(t(sss),breaks=brks,col=colorRampPalette(c("red","black","green"))(numbreaks-1),Rowv=F,Colv=F,trace='none',margins=c(15,15),main=paste0(region[r,1],":",region[r,2],"-",region[r,3]))
    #
    # }
  }

}
