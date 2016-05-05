bgSort <- function( bgfiles, threads=getOption("threads",1L), decreasing=FALSE ){

	ext<-file_ext(bgfiles)
	outnames<-paste0(basename(removeext(bgfiles)),"_scoresort.tmp")
	cmdString <- paste0("sort -T . -k4,4n",if(decreasing){"r"}," ",bgfiles," -o ",outnames)
	res <- cmdRun(cmdString,threads=threads)
	return(outnames)

}
