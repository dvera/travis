intervalBox <-
  function( x , y , breaks=seq(0,1,0.1) , ylims=NULL , quantiles=TRUE , xlabel=NULL , ylabel=NULL , cores="max" , xaxislabel=TRUE , yaxislabel=TRUE  ){

	#check for infinite values? or just discard?
	xinfind <- which(is.infinite(x))
	yinfind <- which(is.infinite(y))
	infind <- unique(c(xinfind,yinfind))

	if(length(infind) > 0){
		x<-x[-infind]
		y<-y[-infind]
	}

	if(is.null(xlabel)){xlabel=deparse(substitute(x))}
	if(is.null(ylabel)){ylabel=deparse(substitute(y))}

  #define score ranges
	if(quantiles){breaks<-quantile(y,probs=breaks)}

	cat("binning x-axis data\n")
	xbinind<-cut(y,breaks,labels=F,include.lowest=T)

	cat("binning y-axis data\n")
	h<-lapply(1:length(breaks),function(i){
		x[which(xbinind==i)]
	})

  boxplot(h,outline=F,if(!is.null(ylims) ) {ylim=ylims})

	binmat<-data.matrix(as.data.frame(mclapply(h,"[[",2)))
	#binmat<-binmat[nrow(binmat):1,]
	row.names(binmat)<-xbreaks[1:breaks]
	colnames(binmat)<-ybreaks[1:breaks]
	binmat<-t(binmat)

	if(white.background){ binmat[binmat==0]<-NA }

	cat("plotting heatmap\n")
	if(grepl("%",maxscore)){maxscore<-quantile(binmat,probs=as.numeric(gsub("%","",maxscore))/100,na.rm=TRUE)}
	if(maxscore==0){maxscore=max(binmat)}
	brks<-seq(0,maxscore,maxscore/100)
	binmat[binmat>brks[101]]<-brks[101]

	if(basic==F){
		layout(matrix(c(4,2,1,3),nrow=2),heights=c(1,4),widths=c(4,1))
		par(oma=c(2,2,2,2))
		par(mar=c(3,3,0,0))
		image(matrix(brks),col=rb(length(brks)-1),breaks=brks,axes=F)
		axis(side=1,at=c(0,1),labels=c(brks[1],brks[length(brks)]))
		mtext("bin counts",side=1,line=2)
		par(mar=c(3,3,0,0))
	}

	image(data.matrix(binmat),breaks=brks,col=rb(length(brks)-1),xaxt='n',yaxt='n')
	text(0,0.9,paste("r =",round(cor(x,y,method="spearman"),3)),pos=4,col="white")
	if(xaxislabel){
		axis(side=1,at=(0:4)/4,labels=xlims[1]+0:4*((xlims[2]-xlims[1])/4))
		mtext(xlabel,side=1,line=2)

	}

	if(yaxislabel){
		axis(side=2,at=(0:4)/4,labels=ylims[1]+0:4*((ylims[2]-ylims[1])/4))
		mtext(ylabel,side=2,line=2)

	}




	if(basic==F){
		dx<-density(x,from=xlims[1],to=xlims[2],na.rm=TRUE)
		dy<-density(y,from=ylims[1],to=ylims[2],na.rm=TRUE)
		par(mar=c(3,0,0,0))
		plot(dy[["y"]],dy[["x"]],type="l",yaxs='i',axes=F,lwd=4)
		par(mar=c(0,3,1,0))
		plot(dx[["x"]],dx[["y"]],type="l",xaxs='i',axes=F,lwd=4,main="")
	}



}
