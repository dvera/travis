alleleSpecificHic <- function( fastqFiles1 , fastqFiles2=NULL , index1prefix, index2prefix , minAS="AS:i:-20", minQual=20 , threads=getOption("threads",1L) , sortBuffer="1G" , sortThreads=NULL , ... ){


  sam_g1r1 <- bowtie2(fastqFiles1, index1prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
  sam_g1r2 <- bowtie2(fastqFiles2, index1prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
  sam_g2r1 <- bowtie2(fastqFiles1, index2prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
  sam_g2r2 <- bowtie2(fastqFiles2, index2prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )


  fields=c("AS","XM","XO","XG","NM")

  pre1 <- paste0(removeext(sam_g1r1),"_preparsed.sam")
  pre2 <- paste0(removeext(sam_g1r2),"_preparsed.sam")

  sys1 <- paste0(removeext(sam_g1r1),"_preparsed.sh")
  sys2 <- paste0(removeext(sam_g1r2),"_preparsed.sh")




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
  PS1=4
  PS2=4+nf

  tag0="\"PG:i:0\""
  tag1="\"PG:i:1\""
  tag2="\"PG:i:2\""
  tag3="\"PG:i:3\""

  fltr <- paste0("for(i=12;i<=NF;i++) { if($i ~ \"",fields,":\"){",fields,"=i } }; if(!",fields,"){",fields,"=\".\"}",collapse=";")

  squareString <- paste(
    "awk -F'\\t' '{",
      "if($1 ~ \"@\"){",
        "print $0}",
      "else if($3==\"*\"){",
        "print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0(rep("\".\"",length(fields)),collapse=","),"}",
      "else{",
        fltr,"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0("$",fields,collapse=","),"} ;",
      paste0(fields,"=\"\"",collapse=";"),
    "}' OFS='\\t' "
  )

  cmdString <- paste0(

      "paste <(",squareString,c(sam_g1r1,sam_g1r2),") <(",squareString,c(sam_g2r1,sam_g2r2),")",
      " | awk -F'\\t' '{",
        "if($1 ~ \"@\"){",
          "print $0 > \"",c(pre1,pre2),"\"}",
        "else if($1!=$",1+nf,"){",
          "print \"ERROR: SAMS DONT SEEM TO BE CONSISTENT\"; exit 1} ",
        "else if($",CH2,"==\"*\" && $",AS1," >= \"",minAS,"\" && $",MQ1," >= ",minQual,"){",
          "print     ",paste0("$",1:nf,collapse=","),",",tag1," > \"",c(pre1,pre2),"\"}",
        "else if($",CH1,"==\"*\" && $",AS2," >= \"",minAS,"\" && $",MQ2," >= ",minQual,"){",
          "print ", paste0("$",nf+(1:(nf)),collapse=","),",",tag2," > \"",c(pre1,pre2),"\"}",
        "else if($",NM1," < $",NM2," && $",XM1," < $",XM2," && $",AS1," >= \"",minAS,"\" && $",MQ1," >= ",minQual," && $",MQ2," >= ",minQual,"){",
          "print    ", paste0("$",1:nf,collapse=","),",",tag1," > \"",c(pre1,pre2),"\"}",
        "else if($",NM2," < $",NM1," && $",XM2," < $",XM1," && $",AS2," >= \"",minAS,"\" && $",MQ1," >= ",minQual," && $",MQ2," >= ",minQual,"){",
          "print ", paste0("$",nf+(1:(nf)),collapse=","),",",tag2," > \"",c(pre1,pre2),"\"}",
        "else if($",MQ1," >= ",minQual," && $",MQ2," >= ",minQual," && $",CH1,"==$",CH2," && ($",PS1,"-$",PS2,")^2 < 100000){",
            "print ", paste0("$",nf+(1:(nf)),collapse=","),",",tag3," > \"",c(pre1,pre2),"\"}",
        "else{",
          "print    ",paste0("$",1:nf,          collapse=","),",",tag0," > \"",c(pre1,pre2),"\"",
        "}",
      "}' OFS='\\t'"
  )
  numfiles <- length(fastqFiles1)
  dump <- lapply(1:numfiles,function(x) write(cmdString[x],sys1[x]))
  dump <- lapply(1:numfiles,function(x) write(cmdString[numfiles+x],sys2[x]))

  cmds <- paste("bash",c(sys1,sys2))
  threads2=floor(threads/3)
  if(threads2<1){threads2<-1}
  cmdRun(cmds,threads=threads2)

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
