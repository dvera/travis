bed12ToBedpe <- function( bedpeFiles , threads=getOption("threads",1L) ){

  outnames <- paste0(basename(removeext(bedpeFiles)),".bed12")

  cmdString <- paste(
    "awk '{if($1==$4 && $2<$5){",
      "print $1,$2,$3,$7,0,\"+\",$2,$6,$8,2,$3-$2\",\"$6-$5,0\",\"$5-$2",
    "} else if ($1==$4 && $2>$5){",
      "print $1,$5,$6,$7,0,\"+\",$6,$2,$8,2,$6-$5\",\"$3-$2,0\",\"$2-$5",
    "}}' OFS='\t' ",bedpeFiles, "> ", outnames
  )

  res <- cmdRun(cmdString,threads=threads)

  return(outnames)
}
