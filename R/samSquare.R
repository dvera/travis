samSquare <- function( samFiles , fields=c("AS","XM","XO","XG","NM") , printHeader=TRUE , sortByName=FALSE , sortBuffer="1G" , sortThreads=NULL , threads=getOption("threads",1L) ){




  if(sortByName){
    outnames <- paste0(basename(removeext(samFiles)),"_square_nsort.sam")
  } else{
    outnames <- paste0(basename(removeext(samFiles)),"_square.sam")
  }

  fltr <- paste0("for(i=12;i<=NF;i++) { if($i ~ \"",fields,":\"){",fields,"=i } }",collapse=";")

  cmdString <- paste(
    "awk -F'\t' '{ if($1 ~ \"@\"){",
      if(printHeader){"print $0"},"}",
      "else if($3==\"*\"){",fltr,"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0(rep("\".\"",length(fields)),collapse=","),"}",
      "else{",fltr,"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0("$",fields,collapse=","),"}",
    "}' OFS='\t'",samFiles,
      if(sortByName){paste("| sort -T . -k1,1 -S",sortBuffer,if(!is.null(sortThreads)){ paste0( "--parallel=" , sortThreads )})} ,
      ">",outnames)

  res <- cmdRun(cmdString, threads=threads)

  return(outnames)

}
