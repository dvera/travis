#' Sort sam files by query name in lexicographic order
#'
#' \code{samSortByName} sorts files by query name in lexicographic order.
#'
#' @param samFiles A character vector of paths to sam files.
#' @param sortThreads A positive integer specifying the number of sorting and compression threads
#' @param threads A positive integer specifying how many bams to process simultaneously.
#' @param memory String specifying maximum memory per thread; suffix K/M/G recognized.
#' @param sortByName Boolean. If TRUE, reads are sorted by name. If FALSE, reads are sorted by chromosome/position


samSortByName <- function( samFiles , threads=getOption("threads",1L), sortBuffer="1G" , sortThreads=NULL ){
	ext<-file_ext(samFiles)
  outnames <- paste0(basename(removeext(samFiles)),"_nsort.sam")
  cmdString <- paste("bash -c 'cat <(grep \"^@\"",samFiles,") <( grep -v \"^@\"",samFiles,"| sort -T . -k1,1 -S",sortBuffer,if(!is.null(sortThreads)){ paste0( "--parallel=" , sortThreads )},") >",outnames,"'")
	res <- cmdRun(cmdString,threads)
	return(outnames)
}
