bgStats <- function( bgFiles , threads=getOption("threads",1L) ){

  cmdString <- paste(
    "awk '{sum+=$4} END{print sum; print NR}' OFS='\t'",bgFiles
  )

  d <- as.data.frame(t(as.data.frame(lapply(cmdRun(cmdString,threads,lines=T),as.numeric))))
  colnames(d)<-c("sum","rows")
  rownames(d)<-basename(bgFiles)
  d$mean <- d$sum/d$rows

  cmdString <- paste(
    "awk '{sumsqdiff+=($4-",d$mean,")^2} END{print sqrt(sumsqdiff/NR)}'", bgFiles
  )

  d$stdev <- unlist( cmdRun(cmdString,threads,intern=T) )

  return(d)
  
}
