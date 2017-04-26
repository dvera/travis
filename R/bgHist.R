#' Plot a histogram of scores sizes in bedGraph files
#'
#' \code{bgHist} plots histograms, density functions, or CDFs of scores in bedGraph files.
#'
#' @param bgFiles A character vector of paths to bedgraph files.
#' @param threads A positive integer specifying how many bedGraph files to process simultaneously.

bgHist <-
function( bedFiles , sample=NULL , chrom = NULL , first = NULL , threads=getOption("threads",1L), ... ){
	options(scipen=9999)
	scores <- bgScores(bedFiles, threads=threads, sample=sample, chrom=chrom, first=first)
	rageHist(scores, threads=threads, ... )
}
