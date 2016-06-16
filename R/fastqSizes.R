fastqSizes <- function( fastqFiles , threads=getOption("threads",1L) ){

  exts <- file_ext(fastqFiles)
  stopifnot(length(unique(exts))==1)

  cmdString <- paste(
    if(exts[1]=="gz"){"zcat"} else{"cat"},
    fastqFiles,
    "| sed 's/\\s/_/g'",
    "| paste - - - -",
    "| cut -f 2",
    "| awk '{print length($0)}'"
  )

  res <- lapply(cmdRun(cmdString,lines=TRUE),as.numeric)
  names(res) <- basename(removeext(fastqFiles))


}
