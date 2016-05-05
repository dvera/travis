matgrid <-
function(mat){
	png(file="test.png",height=5000,width=1000)
	m<-matrix(1:3,nrow=3)
	mat[which(mat>5)]<-5
	mat[which(mat<1)]<-1
	rb<-colorRampPalette(c("black","red"))
	layout(m,heights=c(1,20,5),widths=1)
	par(oma=c(10,0,10,0))
	par(mar=c(10,15,5,15))
	image(matrix(1:100),col=rb(100),breaks=0:100,axes=F,bty="o")
	axis(side=1,at=c(0,1),labels=c(0,50),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
	mtext("sample1_0-200_scg3",side=3,cex=8,outer=T)
	par(mar=c(5,15,5,15),xpd=TRUE)
	image(t(mat),breaks=(1:50)/10,col=rb(49),axes=FALSE)
	lines(rep(-0.05,2),c(0,0.5),lwd=70,col="blue",lend=1)
	axis(side=1,at=c(0,0.5,1),labels=F,lwd=10,padj=1,line=1,tcl=-3)
	par(mar=c(20,15,0,15))
	plot(1:200,colMeans(mat,na.rm=T),type="l",lwd=10,xaxs="i",axes=FALSE)
	axis(side=1,at=c(1,100,200),labels=c(-1000,0,1000),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
	axis(side=2,cex.axis=8,lwd=10,tcl=-3,padj=-1)
	mtext("Distance from TSS (bp)",side=1,cex=5,outer=T,line=2)
	dev.off()
}
