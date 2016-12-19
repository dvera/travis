bed12ToBedpe <- function( bedpeFiles , threads=getOption("threads",1L) ){


  outnames <- paste0(basename(removeext(bedpeFiles)),".bedpe")

  if(any(file.exists(outnames))){stop("file exists")}

  cmdString <- paste(
    "awk '{split($11,sizes,\",\");split($12,starts,\",\");",
      "print $1,$2,$2+sizes[1],$1,$2+starts[2],$2+starts[2]+sizes[2],$4,$5,$6,$6",
    "}' OFS='\t' ",bedpeFiles, "> ", outnames
  )

  res <- cmdRun(cmdString,threads=threads)

  return(outnames)
}
