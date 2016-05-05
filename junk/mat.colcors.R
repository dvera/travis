mat.colcors <-
function( matfile1,matfile2,center="TSS",continue=FALSE,color="black"){
	
	mat1<-data.matrix(read.table(paste(matfile1),stringsAsFactors=FALSE,row.names=NULL))
	mat2<-data.matrix(read.table(paste(matfile2),stringsAsFactors=FALSE,row.names=NULL))
	correlations<-as.numeric(vector(length=ncol(mat1)))
	for(i in 1:ncol(mat1)){
		correlations[i]<-cor(mat1[,i],mat2[,i],use="p")
	}

	if(continue==FALSE){plot(((1:ncol(mat1))-(ncol(mat1)/2)),correlations,xlab=paste("Distance from",center),ylab="Correlation",type="l",col=color,main=paste("Correlation between",removeext(matfile1),"and",removeext(matfile2),"along",center),ylim=c(0,1))}
	else{lines(((1:ncol(mat1))-(ncol(mat1)/2)),correlations,xlab=paste("Distance from",center),ylab="Correlation",col=color,main=paste("Correlation between",removeext(matfile1),"and",removeext(matfile2),"along",center),ylim=c(0,1))}
}
