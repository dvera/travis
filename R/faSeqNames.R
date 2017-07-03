faRemoveChrom <- function( fastaFiles, threads=getOption("threads",1L) ){

  cmdString <- paste0("grep \"^>\" ",fastaFiles," | sed 's/^>//g'")
  res <- cmdRun(cmdString,lines=T,threads=threads)
  
}
