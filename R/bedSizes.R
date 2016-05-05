#' Retreive interval sizes in bed files
#'
#' \code{bedSizes} returns a list of vectors containing the sizes for intervals in a set of bed files.
#'
#' @param bedFiles A character vector of paths to bed files.
#' @param threads A positive integer specifying how many bams to process simultaneously.

bedSizes <- function( bedFiles , sample=NULL , threads=getOption("threads", 1L) ){
  options(scipen=99999)
  numbeds <- length(bedFiles)
  if(is.null(sample)){
    cmdString <- paste("awk '{print $3-$2}'",bedFiles)
  } else{
    stopifnot(is.numeric(sample), length(sample)==1)
    cmdString <- paste("shuf -n",sample,bedFiles,"| awk '{print $3-$2}'")
  }
  res <- cmdRun( cmdString, threads, lines=TRUE)
  res <- lapply( res, as.numeric )
  names(res) <- basename(bedFiles)

  return(res)

}
