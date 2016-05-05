fastq.ssr <- function( fileName , repeatString , minRepeats , outname=NA, outputFasta=TRUE ){

  library(gtools)
  ext<-file_ext(fileName)
  if(is.na(outname)){outname <- paste0(basename(removeext(fileName)), "_", repeatString, minRepeats, ".", if(outputFasta){"fa"} else{ext} ) }

  reg <- unlist(strsplit(repeatString,""))
  reg <- paste0("[",reg)
  reg <- paste0(reg,"]")
  reg <- paste0(reg,collapse="")

  cmdString <- paste(
    "cat",fileName,"| tr '\\t' ' '",
    "| paste - - - - | sed -n",
    paste0("'/\\(",reg,"\\)\\1\\{",minRepeats,",\\}/p'"),
    if(outputFasta){"| cut -f 1,2 | sed 's/@/>/g'"},
    "| tr '\\t' '\\n'",
    ">",outname
  )

  print(cmdString)
  system(cmdString)
  return(outname)
}

fasta.ssr.splice <- function( fileName , repeatString , minRepeats , outname=NA ){

  if(is.na(outname)){outname <- paste0(basename(removeext(fileName)), "_splice", repeatString, minRepeats, ".fa" )}

  reg <- unlist(strsplit(repeatString,""))
  reg <- paste0("[",reg)
  reg <- paste0(reg,"]")
  reg <- paste0(reg,collapse="")

  cmdString <- paste(
    "cat",fileName,
    #"| tr '\\t' ' '",
    "| paste - - | sed",
    paste0("'s/\\(",reg,"\\)\\1\\{",minRepeats,",\\}//g'"),
    "| tr '\\t' '\\n'",
    ">",outname
  )

  print(cmdString)
  system(cmdString)
  return(outname)
}
