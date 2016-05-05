bedRemoveChrom <- function( bedFiles , chroms , outnames=NULL , threads=getOption("threads",1L) ){

  ext <- file_ext(bedFiles)

  if(is.null(outnames)){
    outnames <- paste0(basename(removeext(bedFiles)),"_rmChrom.",ext)
  }

  cmdString <- paste0(
    "grep -v \"",paste0("^",chroms,"[[:space:]]",collapse="\\|"),"\" ",bedFiles," > ",outnames
  )
  res <- cmdRun(cmdString,threads)

  return(outnames)

}
