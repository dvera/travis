bg.discardOverlaps <- function( bgFiles , threads=getOption("threads",1L) ){

  outnames <- paste0( basename(removeext(bgFiles)) , "_bgfix.bg" )

  cmdString <- paste(
    "awk 'BEGIN{a=0}; { if( $2 >= a ){print $0} ; a=$3}' OFS='\t'",bgFiles,">",outnames
  )

  res <- rage.run( cmdString, threads)
  return(outnames)
}
