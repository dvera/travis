#' Sort beds by chromosome and position
#'
#' \code{bedSort} sorts bed files chromosome then position using the system's \code{sort}.
#'
#' @param bedFiles A character vector of paths to bed files.
#' @param sortThreads A positive integer specifying the number of sorting and compression threads
#' @param threads A positive integer specifying how many bams to process simultaneously.
#' @param memory String specifying maximum memory per thread; suffix K/M/G recognized.
#' @param sortByName Boolean. If TRUE, reads are sorted by name. If FALSE, reads are sorted by chromosome/position


bedSort <- function( bedFiles , threads=getOption("threads",1L), sortBuffer="1G" , sortThreads=NULL ){
	ext<-file_ext(bedFiles)
	cmdString <- paste("sort -T . -k1,1 -k2,2n -S",sortBuffer,if(!is.null(sortThreads)){ paste0( "--parallel=" , sortThreads )},bedFiles,"-o",bedFiles)
	res <- cmdRun(cmdString,threads)
	return(bedFiles)
}
