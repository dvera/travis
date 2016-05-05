findcontigs <- function( bedfile, mergewithin = 50, minsize = 500 ){
  mergebed <- bed.merge(bedfile, flank = mergewithin )
  mergebed <- bed.parselengths( mergebed, brks=c( minsize, Inf ) )
  return( mergebed )
}
