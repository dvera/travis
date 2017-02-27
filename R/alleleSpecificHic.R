alleleSpecificHic <- function( fastqFiles1 , fastqFiles2=NULL , index1prefix, index2prefix , minAS="AS:i:-20", minQual=20 , threads=getOption("threads",1L) , sortBuffer="1G" , sortThreads=NULL , ... ){


  sam_g1r1 <- bowtie2(fastqFiles1, index1prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
  sam_g1r2 <- bowtie2(fastqFiles2, index1prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
  sam_g2r1 <- bowtie2(fastqFiles1, index2prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
  sam_g2r2 <- bowtie2(fastqFiles2, index2prefix, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )

  fname1 <- basename(removeext(fastqFiles1))
  fname2 <- basename(removeext(fastqFiles2))

  fields=c("AS","XM","XO","XG","NM")

  i1p <- basename(index1prefix)
  i2p <- basename(index2prefix)

  pre1 <- paste0(fname1,"_genflagged.sam")
  pre2 <- paste0(fname2,"_genflagged.sam")

  sys1 <- paste0(fname1,"_parse.sh")
  sys2 <- paste0(fname2,"_parse.sh")


  par11 <- paste0(fname1,"_",i1p,"-",i1p,".sam")
  par22 <- paste0(fname1,"_",i2p,"-",i2p,".sam")
  par12 <- paste0(fname1,"_",i1p,"-",i2p,".sam")
  par13 <- paste0(fname1,"_",i1p,"-amb.sam")
  par23 <- paste0(fname1,"_",i2p,"-amb.sam")
  par01 <- paste0(fname1,"_",i1p,"-unm.sam")
  par02 <- paste0(fname1,"_",i2p,"-unm.sam")
  par00 <- paste0(fname1,"_unm-unm.sam")
  par33 <- paste0(fname1,"_amb-amb.sam")



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

  # res <- cmdRun(cmdString,threads=threads)

  cmdString <- paste0(
    "paste ",pre1," ",pre2," | awk '{",
      "if($1~\"@\"){",
        paste("print $1,$2,$3 > \"",par11,"\";print $1,$2,$3 > \"",par22,"\";print $1,$2,$3 > \"",par33,"\";print $1,$2,$3 > \"",par00,"\";print $1,$2,$3 > \"",par12,"\";print $1,$2,$3 > \"",par23,"\";print $1,$2,$3 > \"",par01,"\";print $1,$2,$3 > \"",par02,"\";print $1,$2,$3 > \"",par13,"\""),"}",
      "else if($",nf+1,"==\"PG:i:1\" && $",(nf*2)+2,"==\"PG:i:1\"){",
        "P11++; print ",paste0("$",1:nf,collapse=",")," > \"",par11,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par11,"\"}",
      "else if($",nf+1,"==\"PG:i:2\" && $",(nf*2)+2,"==\"PG:i:2\"){",
          "P22++; print ",paste0("$",1:nf,collapse=",")," > \"",par22,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par22,"\"}",
      "else if($",nf+1,"==\"PG:i:3\" && $",(nf*2)+2,"==\"PG:i:3\"){",
          "P33++; print ",paste0("$",1:nf,collapse=",")," > \"",par33,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par33,"\"}",
      "else if($",nf+1,"==\"PG:i:0\" && $",(nf*2)+2,"==\"PG:i:0\"){",
        "P00++; print ",paste0("$",1:nf,collapse=",")," > \"",par00,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par00,"\"}",
      "else if( ($",nf+1,"==\"PG:i:1\" && $",(nf*2)+2,"==\"PG:i:3\") || ($",nf+1,"==\"PG:i:3\" && $",(nf*2)+2,"==\"PG:i:1\")){",
        "P13++; print ",paste0("$",1:nf,collapse=",")," > \"",par13,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par13,"\"}",
      "else if( ($",nf+1,"==\"PG:i:2\" && $",(nf*2)+2,"==\"PG:i:3\") || ($",nf+1,"==\"PG:i:3\" && $",(nf*2)+2,"==\"PG:i:2\")){",
        "P23++; print ",paste0("$",1:nf,collapse=",")," > \"",par23,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par23,"\"}",
      "else if( ($",nf+1,"==\"PG:i:1\" && $",(nf*2)+2,"==\"PG:i:2\") || ($",nf+1,"==\"PG:i:2\" && $",(nf*2)+2,"==\"PG:i:1\")){",
        "P12++; print ",paste0("$",1:nf,collapse=",")," > \"",par12,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par12,"\"}",
      "else if( ($",nf+1,"==\"PG:i:0\" && $",(nf*2)+2,"==\"PG:i:1\") || ($",nf+1,"==\"PG:i:1\" && $",(nf*2)+2,"==\"PG:i:0\")){",
        "P01++; print ",paste0("$",1:nf,collapse=",")," > \"",par01,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par01,"\"}",
      "else if( ($",nf+1,"==\"PG:i:0\" && $",(nf*2)+2,"==\"PG:i:2\") || ($",nf+1,"==\"PG:i:2\" && $",(nf*2)+2,"==\"PG:i:0\")){",
        "P02++; print ",paste0("$",1:nf,collapse=",")," > \"",par02,"\"; print ",paste0("$",(nf+2):((2*nf)+2),collapse=",")," > \"",par02,"\"}",
    "} END{",
      "print \"total mate pairs processed:\",NR;",
      "print \"P11:\",P11,100*P11/NR;",
      "print \"P22:\",P22,100*P22/NR;",
      "print \"P12:\",P12,100*P12/NR;",
      "print \"P13:\",P13,100*P13/NR;",
      "print \"P23:\",P23,100*P23/NR;",
      "print \"P01:\",P01,100*P01/NR;",
      "print \"P02:\",P02,100*P02/NR;",
      "print \"P00:\",P00,100*P00/NR;",
      "print \"P33:\",P33,100*P33/NR",
    "}' OFS='\t'"
  )



}
