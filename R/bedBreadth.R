bedBreadth <- function ( bedfiles ){

  sumcov <- unlist(lapply(bedfiles, function(x) system(paste("awk 'BEGIN{i=0} {i=i+$3-$2} END{print i}'",x),intern=TRUE)))
  return(as.numeric(sumcov))
}
