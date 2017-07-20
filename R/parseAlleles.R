parseAlleles <- function( fastqFiles1 , fastqFiles2=NULL , index1prefix, index2prefix , minAS=-20, minQual=20 , threads=getOption("threads",1L) , fields=c("AS","XM","XO","XG","NM") , ... ){

  #TO DO:
  # pipe bowtie2 to parsing
  # vcf to consensus fastq script
  # output pg tag in final sam files, for PE output


  if(!all(c("AS","XM","XO","XG","NM") %in% fields)){
    stop(paste("the following fields must be present in \"fields\":",paste(c("AS","XM","XO","XG","NM"),collapse=" ") ))
  }

  fname1 <- basename(removeext(fastqFiles1))
  if(index1prefix==index2prefix){
    stop("index1prefix and index2prefix must be different bowtie2 indices\n")
  }

  i1p <- basename(index1prefix)
  i2p <- basename(index2prefix)

  if(i1p==i2p){
    cat("WARNING: index prefix 1 and 2 are identical\n  using g1 for the index found at:",index1prefix,"\n  using g2 for the index found at:",index2prefix,"\n")
  }

  numfiles <- length(fastqFiles1)
  sys1 <- paste0(fname1,"_parse.sh")
  #sys2 <- paste0(fname2,"_parse.sh")

  if(!is.null(fastqFiles2)){
    paired=TRUE
    if(length(fastqFiles1)!=length(fastqFiles2)){
      stop("length of fastqFiles1 and fastqFiles2 must be equal")
    }
    fname2 <- basename(removeext(fastqFiles2))
  }

  # check to see if alignments are provided instead of fastq files
  if(any(file_ext(c(fastqFiles1,fastqFiles2)) %in% c("bam","sam"))){
    sam=TRUE
    cat("alignments detected as inputs\nchecking to see if g1 and g2 files were processed properly\n")
    # count the number of alignments in g1 and g2 files, should be identical
    saml=samtoolsView(c(fastqFiles1,fastqFiles2),count=T,threads=threads)
    sam1l=saml[1:numfiles]
    sam2l=saml[(numfiles+1):(numfiles*2)]

    # get the read names of a chunk of the alignment files, should be identical and in the same order
    sam1h=as.data.frame(cmdRun(paste("samtools view",fastqFiles1,"| head -n 1000 | tail | cut -f 1" ),lines=TRUE),stringsAsFactors=FALSE)
    sam2h=as.data.frame(cmdRun(paste("samtools view",fastqFiles2,"| head -n 1000 | tail | cut -f 1" ),lines=TRUE),stringsAsFactors=FALSE)

    if( !identical(sam1l,sam2l) || !identical(sam1h,sam2h) ){
      stop("g1 and g2 alignments were not processed properly, use fastq files as input")
    } else{
      sam1 <- fastqFiles1
      sam2 <- fastqFiles2
    }

    if(length(unique(sam1h[,1]))<length(sam1h[,1])){
      paired=TRUE
    } else{
      paired=FALSE
    }

  } else{
    sam=FALSE
  }


  if(!sam){
    if(!is.null(fastqFiles2)){
      paired=TRUE
    }
    sam1 <- bowtie2(fastqFiles1, index1prefix, if(paired){fastqFiles2}, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
    sam2 <- bowtie2(fastqFiles1, index2prefix, if(paired){fastqFiles2}, discordant=TRUE, appendIndexToName=TRUE, reorder=TRUE, threads=threads )
  }


  par0 <- paste0("\"",fname1,"_unmapped.sam\"")
  par1 <- paste0("\"",fname1,"_parsed_",i1p,".sam\"")
  par2 <- paste0("\"",fname1,"_parsed_",i2p,".sam\"")
  par3 <- paste0("\"",fname1,"_ambiguous.sam\"")

  par11 <- paste0("\"",fname1,"_",i1p,"-",i1p,".sam\"")
  par22 <- paste0("\"",fname1,"_",i2p,"-",i2p,".sam\"")
  par12 <- paste0("\"",fname1,"_",i1p,"-",i2p,".sam\"")
  par13 <- paste0("\"",fname1,"_",i1p,"-amb.sam\"")
  par23 <- paste0("\"",fname1,"_",i2p,"-amb.sam\"")
  par01 <- paste0("\"",fname1,"_",i1p,"-unm.sam\"")
  par02 <- paste0("\"",fname1,"_",i2p,"-unm.sam\"")
  par00 <- paste0("\"",fname1,"_unm-unm.sam\"")
  par33 <- paste0("\"",fname1,"_amb-amb.sam\"")


  minASstring <- paste0("AS:i:",minAS)
  nf=11+length(fields)
  XM1=paste0("$",which(fields=="XM")+11    )
  XM2=paste0("$",which(fields=="XM")+11+nf )
  AS1=paste0("$",which(fields=="AS")+11    )
  AS2=paste0("$",which(fields=="AS")+11+nf )
  NM1=paste0("$",which(fields=="NM")+11    )
  NM2=paste0("$",which(fields=="NM")+11+nf )
  CH1="$3"
  CH2=paste0("$",3+nf)
  MQ1="$5"
  MQ2=paste0("$",5+nf)
  PS1="$4"
  PS2=paste0("$",4+nf)

  tag0="\"pg:i:0\""
  tag1="\"pg:i:1\""
  tag2="\"pg:i:2\""
  tag3="\"pg:i:3\""

  r1fields  <- paste0("$",1:nf,collapse=",")
  r2fields  <- paste0("$",(nf+2):((2*nf)+2),collapse=",")
  r2fieldsa <- paste0("$",nf+(1:(nf)),collapse=",")
  pgfield2  <- paste0("$",(nf*2)+2)
  pgfield1  <- paste0("$",nf+1)
  minfields1 <- paste0("$",1:11)
  minfields2 <- paste0("$",11+nf+(1:11))



  # this command removes all but the necessary optional fields for parsing alleles
  squareString <- paste(
    "  awk -F'\\t' '{\n",
    "    if($1 ~ \"@\"){\n",
    "      gsub(\"\\t\",\"___\",$0)\n",
    "      print $0\n",
    "    } else if($3==\"*\"){\n",
    "      print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0(rep("\".\"",length(fields)),collapse=","),"\n",
    "    } else{\n",
      paste0(
    "      for(i=12;i<=NF;i++) {\n",
    "        if($i ~ \"",fields,":\"){\n",
    "          ",fields,"=i\n",
    "        }\n",
    "      }\n",
    "      if(!",fields,"){\n",
    "        ",fields,"=\".\"\n",
    "      }\n",
      collapse=""),
    "      print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,",paste0("$",fields,collapse=","),"\n",
    "    }\n",
    "    ",paste0(fields,"=\"\"",collapse=";"),"\n",
    "  }' OFS='\\t' "
  )

  cmdString <- paste0(

        "paste <(\n",
        squareString,sam1,"\n",
        ") <(\n",
        squareString,sam2,") | awk -F'\\t' '{\n",
        "  if($1 ~ \"@\"){\n",
        "    print $1\n",
        "  } else if($1!=$",1+nf,"){\n",
        "    print \"ERROR: SAMS DONT SEEM TO BE CONSISTENT\"\n",
        "    exit 1\n",
        # if g1 and g2 dont align at all
        "  } else if(",CH2,"==\"*\" && ",CH1,"==\"*\"){\n",
        "    print     ",r1fields,",",tag0,"\n",
        # if g1 aligns well but g2 doesnt at all
        "  } else if(",CH2,"==\"*\" && ",AS1," >= \"",minASstring,"\" && ",MQ1," >= ",minQual,"){\n",
        "    print     ",r1fields,",",tag1,"\n",
        # if g2 aligns well but g1 doesnt at all
        "  } else if(",CH1,"==\"*\" && ",AS2," >= \"",minASstring,"\" && ",MQ2," >= ",minQual,"){\n",
        "    print ", r2fieldsa,",",tag2,"\n",
        # if g1 and g2 align well, but the g1 alignment is better
        "  } else if(",NM1," < ",NM2," && ",XM1," < ",XM2," && ",AS1," >= \"",minASstring,"\" && ",MQ1," >= ",minQual," && ",MQ2," >= ",minQual,"){\n",
        "    print    ", r1fields,",",tag1,"\n",
        # if g1 and g2 align well, but the g2 alignment is better
        "  } else if(",NM2," < ",NM1," && ",XM2," < ",XM1," && ",AS2," >= \"",minASstring,"\" && ",MQ1," >= ",minQual," && ",MQ2," >= ",minQual,"){\n",
        "    print ", r2fieldsa,",",tag2,"\n",
        # if g1 and g2 align equally well and to about the same position
        "  } else if(",MQ1," >= ",minQual," && ",MQ2," >= ",minQual," && ",CH1,"==",CH2," && (",PS1,"-",PS2,")^2 < 100000){\n",
        "    print ", r2fieldsa,",",tag3,"\n",
        # if none of the above (e.g. if g1 doesnt align and g2 aligns poorly and/or ambiguously, or g1 and g2 align well but at different positions)
        "  } else{\n",
        "    print ",paste0("$",1:nf, collapse=","),",",tag0,
        "  }\n",
        "}' OFS='\\t' | ",
      # pipe the tagged sam file to the parsing script
        if(paired){
          paste(
        "paste - - | awk -F'\\t' '{\n",
        # if the line is a header
        "    if($1~\"@\"){\n",
        "      HL++\n",
        "      a=gensub(/___/,\"\\t\",g,$1)\n",
        "      print a >",par11,"\n",
        "      print a >",par22,"\n",
        "      print a >",par33,"\n",
        "      print a >",par00,"\n",
        "      print a >",par12,"\n",
        "      print a >",par23,"\n",
        "      print a >",par01,"\n",
        "      print a >",par02,"\n",
        "      print a >",par13,"\n",
            # if r1 and r2 are g1
        "    } else if(",pgfield1,"==",tag1," && ",pgfield2,"==",tag1,"){\n",
        "      P11++\n",
        "      print ",r1fields," >",par11,"\n",
        "      print ",r2fields," >",par11,"\n",
            # if r1 and r2 are g2
        "    } else if(",pgfield1,"==",tag2," && ",pgfield2,"==",tag2,"){\n",
        "      P22++\n",
        "      print ",r1fields," >",par22,"\n",
        "      print ",r2fields," >",par22,"\n",
            # if r1 and r2 are ambiguous
        "    } else if(",pgfield1,"==",tag3," && ",pgfield2,"==",tag3,"){\n",
        "      P33++\n",
        "      print ",r1fields," >",par33,"\n",
        "      print ",r2fields," >",par33,"\n",
            # if r1 and r2 are unmapped
        "    } else if(",pgfield1,"==",tag0," && ",pgfield2,"==",tag0,"){\n",
        "      P00++\n",
        "      print ",r1fields," >",par00,"\n",
        "      print ",r2fields," >",par00,"\n",
            # if r1 is g1 and r2 is ambiguous
        "    } else if( (",pgfield1,"==",tag1," && ",pgfield2,"==",tag3,") || (",pgfield1,"==",tag3," && ",pgfield2,"==",tag1,")){\n",
        "      P13++\n",
        "      print ",r1fields," >",par13,"\n",
        "      print ",r2fields," >",par13,"\n",
            # if r1 is g2 and r2 is ambiguous
        "    } else if( (",pgfield1,"==",tag2," && ",pgfield2,"==",tag3,") || (",pgfield1,"==",tag3," && ",pgfield2,"==",tag2,")){\n",
        "      P23++\n",
        "      print ",r1fields," >",par23,"\n",
        "      print ",r2fields," >",par23,"\n",
            # if r1 is g1 and r2 is g2
        "    } else if( (",pgfield1,"==",tag1," && ",pgfield2,"==",tag2,") || (",pgfield1,"==",tag2," && ",pgfield2,"==",tag1,")){\n",
        "      P12++\n",
        "      print ",r1fields," >",par12,"\n",
        "      print ",r2fields," >",par12,"\n",
            # if r1 is g1 and r2 is unmapped
        "    } else if( (",pgfield1,"==",tag0," && ",pgfield2,"==",tag1,") || (",pgfield1,"==",tag1," && ",pgfield2,"==",tag0,")){\n",
        "      P01++\n",
        "      print ",r1fields," >",par01,"\n",
        "      print ",r2fields," >",par01,"\n",
            # if r1 is g2 and r2 is unmapped
        "    } else if( (",pgfield1,"==",tag0," && ",pgfield2,"==",tag2,") || (",pgfield1,"==",tag2," && ",pgfield2,"==",tag0,")){\n",
        "      P02++\n",
        "      print ",r1fields," >",par02,"\n",
        "      print ",r2fields," >",par02,"\n",
        "    } } END{\n",
        "    print \"total mate pairs processed:\",NR-HL\n",
        "    print \"P11:\",P11,100*P11/(NR-HL)\n",
        "    print \"P22:\",P22,100*P22/(NR-HL)\n",
        "    print \"P12:\",P12,100*P12/(NR-HL)\n",
        "    print \"P13:\",P13,100*P13/(NR-HL)\n",
        "    print \"P23:\",P23,100*P23/(NR-HL)\n",
        "    print \"P01:\",P01,100*P01/(NR-HL)\n",
        "    print \"P02:\",P02,100*P02/(NR-HL)\n",
        "    print \"P00:\",P00,100*P00/(NR-HL)\n",
        "    print \"P33:\",P33,100*P33/(NR-HL)\n",
        "}' OFS='\\t'\n")

      } else{
        paste(
        "awk -F'\\t' '{\n",
          # if the line is a header
        "    if($1~\"@\"){\n",
        "      HL++\n",
        "      $0=$1\n",
        "      a=gensub(/___/,\"\\t\",g,$1)\n",
        "      print a > ",par1,"\n",
        "      print a > ",par2,"\n",
        "      print a > ",par3,"\n",
        "      print a > ",par0,"\n",
        # if read is g1
        "    } else if(",pgfield1,"==",tag1,"){\n",
        "      P1++\n",
        "      print ",r1fields," > ",par1,"\n",
        # if read is g2
        "    } else if(",pgfield1,"==",tag2,"){\n",
        "      P2++\n",
        "      print ",r1fields," > ",par2,"\n",
        # if read is ambiguous
        "    } else if(",pgfield1,"==",tag3,"){\n",
        "      P3++\n",
        "      print ",r1fields," > ",par3,"\n",
        # if read is unmapped
        "    } else if(",pgfield1,"==",tag0,"){\n",
        "      P0++\n",
        "      if($3==\"*\"){\n",
        "        print",minfields1," > ",par0,"\n",
        "      } else{\n",
        "        print ",r1fields," > ",par0,"\n",
        "      }\n",
        "  } } END{\n",
        "    print \"total reads processed:\",NR-HL\n",
        "    print \"P0:\",P0,100*P0/(NR-HL)\"%\"\n",
        "    print \"P1:\",P1,100*P1/(NR-HL)\"%\"\n",
        "    print \"P2:\",P2,100*P2/(NR-HL)\"%\"\n",
        "    print \"P3:\",P3,100*P3/(NR-HL)\"%\" \n",
      "}' OFS='\\t'\n" )
      }

  )

  dump <- lapply(1:numfiles,function(x) write(cmdString[x],sys1[x]))
  cmds <- paste("bash",sys1)
  threads2=floor(threads/4)
  if(threads2<1){threads2<-1}
  cmdRun(cmds,threads=threads2)

}
