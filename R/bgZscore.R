bgZscores <- function( bgFiles , threads=getOption("threads",1L) ){

  outnames <- paste0(basename(removeext(bgFiles)),"_zscore.bg")

  smry <- bgStats(bgFiles, threads=threads)

  cmdString <- paste(
    "awk '{$4=($4-",smry$mean,")/",smry$stdev,";print $0}' OFS='\t'", bgFiles,">",outnames
  )

  res <- cmdRun(cmdString, threads )

  return(outnames)

}
