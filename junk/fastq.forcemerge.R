fastq.forcemerge <- function( r1 , r2 , minLength=50 , minOverlap=10 , pValue=0.01 , memory="10G" , threads="max" ){

  library(parallel)
  if(threads=="max"){threads=detectCores()-1}

  if(length(r1)!=length(r2)){stop("must have equal number of forward and reverse reads")}
  numpairs <- length(r1)

  if(all(grepl("R1",r1))){
    outnames <- gsub("R1","fmerged",basename(r1))
  } else{
    outnames <- paste0(basename(removeext(r1)),"_fmerged.fastq")
  }
  cmdString <- paste(

    "paste",
    r1,
    r2,
    "| awk '{r=NR%4; if(r==2 || r==0){print $1 $2} else{print $1}}' >",
    outnames
  )

  for(i in 1:numpairs){
    print(cmdString[i])
    system(cmdString[i])
  }

  return(outnames)

}
