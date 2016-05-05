track.wig <-
function( trackfiles,tracknames=removeext(trackfiles), plotcolors=rep("black",length(trackfiles)), user=Sys.getenv("USER"),server="epsilon.bio.fsu.edu",path="public_html/hubs/dlv",genome="hg19",range=c(0,50),printtrack=F){
	library(tools)
	if(user=="dlv04c"){user="dvera"}
	if(length(which(file_ext(trackfiles) %in% c("wig","Wig","wiggle")) > 0)){stop("cannot create track from wig files, convert to bigWig")}
	print(data.frame("file"=trackfiles,"names"=tracknames,"color"=plotcolors))
	
	track<-""
	for(i in 1:length(trackfiles)){
		c<-paste(col2rgb(plotcolors[i]),collapse=",")
		subtrack<-c(
			paste("track ",tracknames[i],sep=""),
			paste("bigDataUrl bbi/",trackfiles[i],sep=""),
			paste("shortLabel ",tracknames[i],sep=""),
			paste("longLabel ",tracknames[i],sep=""),
			paste("type bigWig ",range[1]," ",range[2],sep=""),
			paste("color ",c,sep=""),
			paste("altColor ",c,sep=""),

		""
		)
		track<-append(subtrack,track)
	}
	
	if(printtrack==TRUE){cat(unlist(outfile),sep="\n")}
	
	trackfile<-data.frame("V1"=track,stringsAsFactors=FALSE)
	write.tsv(trackfile,file="tmphubdb.txt")
	filelist<-paste(trackfiles,sep=" ",collapse=" ")
	system(paste("scp ",filelist," ",user,"@",server,":",path,"/",genome,"/bbi/",sep=""))
	system(paste("cat tmphubdb.txt | ssh ",user,"@",server," 'cat >> ",path,"/",genome,"/hubDb.txt'",sep=""))
}
