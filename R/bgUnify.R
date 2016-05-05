bgUnify <- function( bgFiles , filler="NA" , discardUnshared=TRUE){
  outnames <- paste0(basename(removeext(bgFiles)),"_unified.bg")
  if(discardUnshared){filler="NA"}
  numfiles <- length(bgFiles)
  cmdString <- paste(
    "bedtools unionbedg",
    "-filler",filler,
    "-i",paste(bgFiles, collapse=" "),
    if(discardUnshared){"| grep -v 'NA'"},
    "| awk '{",
    paste0("print $1,$2,$3,$",3+(1:numfiles),"> \"",outnames,"\"",collapse=";"),
    "}' OFS='\t'"
  )
  print(cmdString)
  system(cmdString)
  return(outnames)

}
