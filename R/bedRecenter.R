bedRecenter <-
function( bed,regionsize,center=FALSE,strand=TRUE,start=2,stop=3){

	bedname<-basename(removeext(bed))
	windowbed<-read.tsv(bed)
	winbedname<-paste0(bedname,"_recenter.bed")

	if(strand==TRUE){
		negrows<-which(windowbed[,6]=="-")
		windowbed[negrows,2]<-windowbed[negrows,stop]
		windowbed[-negrows,2]<-windowbed[-negrows,start]
	} else{
		windowbed[,2]=windowbed[,start]
	}

	if(center==TRUE){
		windowbed[,2]=round((windowbed[,start]+windowbed[,stop])/2)
	}

	windowbed[,2]<-windowbed[,2]-regionsize/2
	windowbed[,3]<-windowbed[,2]+regionsize
	windowbed<-windowbed[which(windowbed[,2]> 0),]
	#windowbed<-windowbed[order(windowbed$V1,windowbed$V2),]
	write.tsv(windowbed,file=winbedname)
	return(winbedname)

}
