bedpeToIbed <- function( bedpeFiles , threads=getOption("threads",1L) ){


  outnames <- paste0(basename(removeext(bedpeFiles)),".ibed")

  if(any(file.exists(outnames))){stop("file exists")}

  cmdString <- paste(
    "awk",
      "print $1,$2,$3,$4 \":\"$5\"-\"$6\",\"",
    "}' OFS='\t' ",bedpeFiles, "> ", outnames
  )

  res <- cmdRun(cmdString,threads=threads)

  return(outnames)
}
