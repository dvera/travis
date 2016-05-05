vplot.draw <-
function( vmat,hm.max="default",png=FALSE, plotcolors=c("black","green")){
	colramp <- colorRampPalette(plotcolors, space = "rgb")
	library(gplots)
	matname<-basename(removeext(vmat))
	m<-read.mat(vmat)
	if(hm.max=="default"){hm.max=as.numeric(round(quantile(m[m!=0],prob=0.99)))}
	filename<-paste(matname,"_max",round(hm.max),"_vplot.png",sep="")
	if(png==TRUE){
		png(file=filename,width=1000,height=1000)
	}
	heatmap.2(m,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=seq(0,hm.max,by=hm.max/100),labRow=NA,labCol=NA,col=colramp(100))
	if(png==TRUE){
		dev.off()
	}
}
