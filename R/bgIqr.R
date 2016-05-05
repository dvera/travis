bgIqr <-
function( bgs, ... ){

  q <- bgQuantiles( bgs, quantiles=c(0.25,0.75) )
  iqrs <- q[2,]-q[1,]
  #names(iqrs) <- NULL
  return(iqrs)

}
