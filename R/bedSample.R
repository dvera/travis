bedSample <-
function( beds , count , threads="max" ){
	library(tools)
	ext<-file(ext(bed))
	#totallines<-filelines(bed)
	#if(count > totallines ){stop("sample quantity is greater than available lines")}
	outname<-paste0(basename(removeext(bed)),"_random",count,".",ext)
	system(paste("shuf -n",count,beds,">",newname))
}
