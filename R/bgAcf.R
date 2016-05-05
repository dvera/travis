#bgAcf <- function ( bgs , maxlag=9 , maintitle="", returnlag=1, plot=TRUE ){
bgAcf <- function ( bgs , maxlag=10, threads=getOption("threads",1L) ){

  numbgs=length(bgs)
  bgl <- bgRead(bgs,makematrix=TRUE,enforceEquality=TRUE,threads=threads)
  bgl <- as.data.frame(bgl)
  acfmat = t(as.data.frame( lapply( lapply( bgl, acf, lag.max=maxlag, plot=F), "[[", 1 ) ))
  #a=-1*data.matrix(as.data.frame((lapply(a,diff))))

  acfmat=acfmat[,-1]
  colnames(acfmat) <- paste0("lag",1:maxlag)
  #acor=a[returnlag,]

  #a=a[,order(a[1,],decreasing=T)]


  # if(plot){
  #   par(mar=c(15,4,4,4),family="mono")
  #   plot(0,type="n",xlim=c(1,numbgs),ylim=c(0,1),ylab="autocorrelation",xlab="",xaxt="n",main=maintitle)
  #   abline(h=(0:10)/10,col=rgb(0,0,0,20,maxColorValue=255) )
  #   abline(v=1:numbgs,col=c(rgb(255,0,0,40,maxColorValue=255),rgb(0,255,0,40,maxColorValue=255)))
  #   for(i in 1:nrow(a)){lines(1:numbgs,a[i,],col=rainbow(nrow(a))[i])}
  #   for(i in 1:nrow(a)){points(1:numbgs,a[i,],col=rainbow(nrow(a))[i],pch=20)}
  #   #abline(v=1:numbgs,col=rgb(0,0,0,40,maxColorValue=255))
  #   axis(1,at=1:numbgs,labels=colnames(a),las=3,cex.axis=0.7)
  #   legend("topright",title="lag",legend=1:maxlag,col=rainbow(maxlag),lwd=3,cex=0.9)
  # }




  return(acfmat)




}
