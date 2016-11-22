bedSample <-function( beds , count , threads=getOption(threads,1L) ){
	library(tools)
	ext<-file_ext(beds)
	#totallines<-filelines(bed)
	#if(count > totallines ){stop("sample quantity is greater than available lines")}
	outname<-paste0(basename(removeext(beds)),"_random",count,".",ext)
	cmdString=paste("shuf -n",count,beds,">",outname)
	res <- cmdRun(cmdString,threads=threads)
	return(outname)
}
