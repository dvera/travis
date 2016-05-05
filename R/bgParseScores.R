bgParseScores <-
function( bgfile, range=c(1,Inf) ){
	library(tools)
	bgname<-basename(removeext(bgfile))
	ext<-file_ext(bgfile)
	if(range[1]==-Inf){
		outname<-paste(bgname,"_lt",range[2],".",ext,sep="")
		system(paste("awk '$4 <=",range[2],"'",bgfile,">",outname))
	}
	if(range[2]==Inf){
		outname<-paste(bgname,"_gt",range[1],".",ext,sep="")
		system(paste("awk '$4 >=",range[1],"'",bgfile,">",outname))
	}
	if(range[1] != -Inf & range[2] != Inf){
		outname<-paste(bgname,"_gt",range[1],"_lt",range[2],".",ext,sep="")
		system(paste("awk '$4 >=",range[1],"'",bgfile,">",outname))
	}
	return(outname)
}
