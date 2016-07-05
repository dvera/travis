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

samParseGenomes <- function( genome1sams , genome2sams , minAS=-20, minQual=20 , threads=getOption("threads",1L) , sortBuffer="1G" ){

  # f = files("*.fastq")
  # ft = cutadapt( f )

  # s1 = bowtie2( ft , index1 )
  # s2 = bowtie2( ft , index2 )

  # check if list of mate pairs given
  if(is.list(genome1sams) | is.list(genome2sams) ){
    if(length(unique(c(length(genome1sams[[1]]),length(genome1sams[[2]]),length(genome2sams[[1]]),length(genome2sams[[2]]) ))) >1 ){stop("length of mate pairs should be equal")}
    numsamples <- length(genome1sams[[1]])
    paired=TRUE
  } else{
    if(length(unique(c(length(genome1sams),length(genome2sams)))) > 1 ){stop("genome sams are not paired")}
    paired=FALSE
    numsamples <- length(genome1sams)
  }

  allsams <- c(unlist(genome1sams),unlist(genome2sams))

  if(paired){
    g1 = 1:(numsamples*2)
    g2 = (2*numsamples+1):(numsamples*4)
    r1 <- c( (1:numsamples),(numsamples*2)+(1:numsamples) )
    r2 <- c( numsamples+(1:numsamples) , (numsamples*3)+(1:numsamples) )
  } else{
    g1 = 1:numsamples
    g2 = (numsamples+1):(numsamples*2)
    r1 <- 1:(numsamples*2)
  }

  # generate some names for temporary and output files
  basenames <- basename(removeext(unlist(allsams)))
  preparsed <- paste0(basenames, "_preparsed.sam")
  unnparsed <- paste0(basenames, "_unnparsed.sam")
  spcparsed <- paste0(basenames, "_spcparsed.sam")
  genparsed <- paste0(basenames, "_genparsed.sam")

  # check if sams have identical query sequences
  cat("checking if sams have identical queries\n")
  samcounts <- matrix(samtoolsView(allsams,count=TRUE,threads=threads),nrow=numsamples)
  samcountstable <- apply(samcounts,1,table)

  if(class(samcountstable)=="integer"){
    unequal=FALSE
    cat("genome sams appear to have the same set of quieries\n")
  } else{
    unequal=TRUE
    cat("WARNING: genome1 and genome2 alignments do not appear to have the same set of queries. Including unaligned sequences and reporting only one alignment for every read will make allele parsing faster.\n")
  }

  cat("sorting reads and standardizing optional flags\n")
  fields=c("AS","XM","XO","XG","NM")

  allsamsSorted = samSquare ( allsams , fields=fields, sortByName=T, printHeader=F, threads=threads, sortBuffer=sortBuffer )

  nf=11+length(fields)

  XM1=which(fields=="XM")+11
  XM2=XM1+nf-1
  AS1=which(fields=="AS")+11
  AS2=AS1+nf-1
  NM1=which(fields=="NM")+11
  NM2=NM1+nf-1

  cat("parsing reads common to both genome sam files\n")
  cmdString <- paste0(
    "join -t $'\t' -j 1 ",allsamsSorted[g1]," ",allsamsSorted[g2], " | awk -F'\t' '{",
      "split($",AS1,",AS1,\":\"); split($",AS2,",AS2,\":\");",
           "if($",2+nf,"==\"*\" && AS1[3] > ",minAS," && $5 > ",minQual,"){ print ",paste0("$",1:nf,collapse=",")," > \"",preparsed[g1],"\"}",
      "else if($3==\"*\" && AS1[3] > ",minAS," && $",4+nf," > ",minQual,"){ print $1,", paste0("$",nf+(1:(nf-1)),collapse=",")," > \"",preparsed[g2],"\"}",
      "else if($",NM1," < $",NM2," && $",XM1," < $",XM2," && AS1[3] > ",minAS," && $5 >= ",minQual," && $",4+nf," > ",minQual,"){ print    ", paste0("$",1:nf,collapse=",")         ," > \"",preparsed[g1],"\"}",
      "else if($",NM2," < $",NM1," && $",XM2," < $",XM1," && AS2[3] > ",minAS," && $5 >= ",minQual," && $",4+nf," > ",minQual,"){ print $1,", paste0("$",nf+(1:(nf-1)),collapse=",")," > \"",preparsed[g2],"\"}",
      "else{",
        "print    ",paste0("$",1:nf,          collapse=",")," > \"",unnparsed[g1],"\";",
        "print $1,",paste0("$",nf+(1:(nf-1)), collapse=",")," > \"",unnparsed[g2],"\"",
      "}",
    "}' OFS='\t'"
  )
  res <- cmdRun(cmdString,threads=threads)

  if(unequal){

    cat("parsing reads unique to genome 1\n")
    cmdString <- paste0(
      "join -v 1 -t $'\t' -j 1 ",allsamsSorted[g1]," ",allsamsSorted[g2], " | awk -F'\t' '{",
        #"if(NR==1){nf=NF};",fltr1,"; split($AS1,AS1S,\":\");",
        "split($",AS1,",AS1,\":\");",
        "if(AS1[3] > ",minAS," && $5 >= ",minQual,"){ print ",paste0("$",1:nf,collapse=",")," > \"",spcparsed[g1],"\"}",
        "else{print $0 >> \"",unnparsed[g1],"\"}",
      "}' OFS='\t'"
    )
    res <- cmdRun(cmdString,threads=threads)

    cat("parsing reads unique to genome 2\n")
    cmdString <- paste0(
      "join -v 1 -t $'\t' -j 1 ",allsamsSorted[g1]," ",allsamsSorted[g2], " | awk -F'\t' '{",
        #"if(NR==1){nf=NF};",fltr1,"; split($AS1,AS1S,\":\");",
        "split($",AS1,",AS1,\":\");",
        "if(AS1[3] > ",minAS," && $5 >= ",minQual,"){ print ",paste0("$",1:nf,collapse=",")," > \"",spcparsed[g2],"\"}",
        "else{print $0 > \"",unnparsed[g2],"\"}",
      "}' OFS='\t'"
    )
    res <- cmdRun(cmdString,threads=threads)

  }


  cmdString <- paste(
    "awk '{if($1~\"@\"){print $0} else{exit 0}}' OFS='\t'",allsams[r1],">", genparsed[r1]
  )
  res <- cmdRun(cmdString,threads=threads)

  if(unequal | paired){

    # if genome-specific sams don't exist, make them so merging doesn't fail.
    if(unequal){
      noe <- which(!file.exists(spcparsed))
      if(length(noe)>1){file.create(spcparsed[noe])}
    }

    cat("consolidating parsed reads\n")
    cmdString <- paste(
      "sort -k1,1 -m",
      preparsed[r1],
      if(paired){preparsed[r2]},
      if(unequal){spcparsed[r1]},
      if(paired & unequal){spcparsed[r2]},
      "| awk '{if($1!=p){print $0}; p=$1}' OFS='\t' >> ",genparsed[r1]
    )
    res <- cmdRun(cmdString,threads=threads)
  } else{
    cmdString <- paste("cat",preparsed[r1],">>",genparsed[r1])
    res <- cmdRun(cmdString, threads=threads)
  }

  return(genparsed[r1])

}
