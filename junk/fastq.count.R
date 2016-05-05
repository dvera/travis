fastq.count <- function( fastqs , threads=getOption("threads",1L) ){
  return(filelines( fastqs , threads )/4)
}
