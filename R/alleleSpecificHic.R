alleleSpecificHic <- function( fastqFiles1 , fastqFiles2=NULL , index1prefix, index2prefix , minAS="AS:i:-20", minQual=20 , threads=getOption("threads",1L) , sortBuffer="1G" , sortThreads=NULL , ... ){



  # genome1sams <- bowtie2(fastqFiles1, fastqFiles2, index1prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads, ... )
  # genome2sams <- bowtie2(fastqFiles1, fastqFiles2, index2prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads, ... )

  sam_g1r1 <- bowtie2(fastqFiles1, index1prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads, ... )
  sam_g1r2 <- bowtie2(fastqFiles2, index1prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads, ... )
  sam_g2r1 <- bowtie2(fastqFiles1, index2prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads, ... )
  sam_g2r2 <- bowtie2(fastqFiles2, index2prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads, ... )


  fields=c("AS","XM","XO","XG","NM")
  fltr <- paste0("for(i=12;i<=NF;i++) { if($i ~ \"",fields,":\"){",fields,"=i } }",collapse=";")

  squareString <- paste(
    "awk -F'\t' '{",
      "if($1 ~ \"@\"){}",
      "else if($3==\"*\"){",fltr,"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0(rep("\".\"",length(fields)),collapse=","),"}",
      "else{",fltr,"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0("$",fields,collapse=","),"}",
    "}' OFS='\t'"
  )

  # check if list of mate pairs given
  # if(is.list(genome1sams) | is.list(genome2sams) ){
  #   if(length(unique(c(length(genome1sams[[1]]),length(genome1sams[[2]]),length(genome2sams[[1]]),length(genome2sams[[2]]) ))) >1 ){stop("length of mate pairs should be equal")}
  #   numsamples <- length(genome1sams[[1]])
  #   paired=TRUE
  # } else{
  #   if(length(unique(c(length(genome1sams),length(genome2sams)))) > 1 ){stop("genome sams are not paired")}
  #   paired=FALSE
  #   numsamples <- length(genome1sams)
  # }

  # allsams <- c(unlist(genome1sams),unlist(genome2sams))
  #
  # if(paired){
  #   g1 = 1:(numsamples*2)
  #   g2 = (2*numsamples+1):(numsamples*4)
  #   r1 <- c( (1:numsamples),(numsamples*2)+(1:numsamples) )
  #   r2 <- c( numsamples+(1:numsamples) , (numsamples*3)+(1:numsamples) )
  # } else{
  #   g1 = 1:numsamples
  #   g2 = (numsamples+1):(numsamples*2)
  #   r1 <- 1:(numsamples*2)
  # }

  # generate some names for temporary and output files
  # basenames <- basename(removeext(unlist(allsams)))
  #
  #
  #
  # preparsed <- paste0(basenames, "_preparsed.sam")
  # unnparsed <- paste0(basenames, "_unnparsed.sam")
  # spcparsed <- paste0(basenames, "_spcparsed.sam")
  # genparsed <- paste0(basenames, "_genparsed.sam")

  pre1 <- paste0(removeext(sam_g1r1),"_preparsed.sam")
  pre2 <- paste0(removeext(sam_g1r2),"_preparsed.sam")
  # check if sams have identical query sequences
  # cat("checking if sams have identical queries\n")
  # samcounts <- matrix(samtoolsView(allsams,count=TRUE,threads=threads),nrow=numsamples)
  # samcountstable <- apply(samcounts,1,table)
  #
  # if(class(samcountstable)=="integer"){
  #   unequal=FALSE
  #   cat("genome sams appear to have the same set of quieries\n")
  # } else{
  #   unequal=TRUE
  #   cat("WARNING: genome1 and genome2 alignments do not appear to have the same set of queries. Including unaligned sequences and reporting only one alignment for every read will make allele parsing faster.\n")
  # }

  #cat("sorting reads and standardizing optional flags\n")
  #fields=c("AS","XM","XO","XG","NM")

  #allsamsSorted = samSquare ( allsams , fields=fields, sortByName=T, printHeader=F, threads=threads, sortBuffer=sortBuffer, sortThreads=sortThreads )

  nf=11+length(fields)

  XM1=which(fields=="XM")+11
  XM2=XM1+nf

  AS1=which(fields=="AS")+11
  AS2=AS1+nf

  NM1=which(fields=="NM")+11
  NM2=NM1+nf

  CH1=3
  CH2=3+nf

  MQ1=5
  MQ2=5+nf

  # XM_G1R1=which(fields=="XM")+11
  # XM_G1R2=XM1+nf
  # XM_G2R1=XM1+nf+nf
  # XM_G2R2=XM1+nf+nf+nf
  #
  # AS_G1R1=which(fields=="AS")+11
  # AS_G1R2=AS1+nf
  # AS_G2R1=AS1+nf+nf
  # AS_G2R2=AS1+nf+nf+nf
  #
  # NM_G1R1=which(fields=="NM")+11
  # NM_G1R2=NM1+nf
  # NM_G2R1=NM1+nf+nf
  # NM_G2R2=NM1+nf+nf+nf
  #
  # CH_G1R1=3
  # CH_G1R2=3+nf
  # CH_G2R1=3+nf+nf
  # CH_G2R2=3+nf+nf+nf
  #
  # MQ_G1R1=5
  # MQ_G1R2=5+nf
  # MQ_G2R1=5+nf+nf
  # MQ_G2R2=5+nf+nf+nf
  #

  tag0="PG\":\"i\":\"0"
  tag1="PG\":\"i\":\"1"
  tag2="PG\":\"i\":\"2"


  cat("parsing reads common to both genome sam files\n")
  cmdString <- paste0(
    #"bash -c paste <(",squareString,genome1sams,"| paste",rep("-",2*nf),") <(",squareString,genome2sams,"| paste",rep("-",2*nf),") | awk -F'\t' '{",
      #"split($",AS_G1R1,",AS_G1R1,\":\"); split($",AS_G1R2,",AS_G1R1,\":\"); split($",AS_G2R1,",AS_G2R1,\":\"); split($",AS_G2R2,",AS_G2R2,\":\");",

      "bash -c paste <(",squareString,sam_g1r1,") <(",squareString,sam_g2r1,")",
      " | awk -F'\t' '{",
        "if($",CH2,"==\"*\" && $",AS1,"[3] > ",minAS," && $",MQ1," > ",minQual,"){",
          "print     ",paste0("$",1:nf,collapse=","),",",tag1," > \"",pre1,"\"}",
        "else if($",CH1,"==\"*\" && $",AS2,"[3] > ",minAS," && $",MQ2," > ",minQual,"){",
          "print $1,", paste0("$",nf+(1:(nf-1)),collapse=","),",",tag2," > \"",pre1,"\"}",
        "else if($",NM1," < $",NM2," && $",XM1," < $",XM2," && $",AS1,"[3] > ",minAS," && $",MQ1," >= ",minQual," && $",MQ2," >= ",minQual,"){",
          "print    ", paste0("$",1:nf,collapse=","),",",tag1," > \"",pre1,"\"}",
        "else if($",NM2," < $",NM1," && $",XM2," < $",XM1," && AS2[3] > ",minAS," && $5 >= ",minQual," && $",4+nf," > ",minQual,"){",
          "print $1,", paste0("$",nf+(1:(nf-1)),collapse=","),",",tag2," > \"",pre1,"\"}",
        "else{",
          "print    ",paste0("$",1:nf,          collapse=","),",",tag0," > \"",pre1,"\"}",
        "}",
      "}' OFS='\t'"
  )
  res <- cmdRun(cmdString,threads=threads)


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
