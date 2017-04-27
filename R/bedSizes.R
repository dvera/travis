#' Retreive interval sizes in bed files
#'
#' \code{bedSizes} returns a list of vectors containing the sizes for intervals in a set of bed files.
#'
#' @param bedFiles A character vector of paths to bed files.
#' @param threads A positive integer specifying how many bams to process simultaneously.

bedSizes <- function( bedFiles , sample=NULL , chrom=NULL , first=NULL , threads=getOption("threads", 1L) ){
  options(scipen=99999)
  numbeds <- length(bedFiles)
  
  if( !is.null(sample) & !is.null(chrom) ){stop("must set sample or chrom or neither, but not both")}

  if(!is.null(sample)){
    stopifnot(is.numeric(sample), length(sample)==1,sample>0)
    cmdString <- paste("shuf -n",sample,bedFiles,"| awk '{print $3-$2}'")
  } else if( !is.null(chrom)){
    stopifnot(length(chrom)==1)
    cmdString <- paste0("grep -P '^",chrom,"\\t' ",bedFiles," | awk '{print $3-$2}'")
  } else if(!is.null(first) & is.numeric(first)){
    cmdString <- paste("head -n",first,bedFiles," | awk '{print $3-$2}'")
  } else{
    cmdString <- paste("awk '{print $3-$2}'",bedFiles)
  }

  res <- cmdRun( cmdString, threads, lines=TRUE)
  res <- lapply( res, as.numeric )
  names(res) <- basename(bedFiles)

  return(res)

}
