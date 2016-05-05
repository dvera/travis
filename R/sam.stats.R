#' @title Get summary statistics for bam files
#'
#' @description Calculates the number of duplicate, unique-aligning, multi-aligning, unaligned, and total reads for a set of position-sorted bam files.
#'
#' @param bamFiles a character vector of paths to position-sorted bam files
#' @param genomefile a string of the path to the genome file, a tab-delimited file with chromosome names in column 1 and chromosome sizes in column 2
#' @param probes a string of the path to a bed file containing the locations of specific regions of interest to additionally calculate the overlap of reads to these regions
#' @param expandprobes a positive integer that is used to expand the region around intervals in 'probes' in which to calculate overlap of reads
#' @param plot boolean indicating if a summary plot is drawn for calculated statistics
#' @param uniquescore a positive integer used as the threshold for calling an alignment unique
#' @param cores a positive integer specifying the number of bam files to process simultaneously

sam.stats <- function (bamFiles, genomefile, probes=NULL , expandprobes=200 , plot=TRUE, uniquescore=10, cores=11){

  # TO DO
  #   check if bam files are sorted
  #   alignment quality cdf

  numfiles=length(bamFiles)


#  cat("\n###########################\n")
  cat("### counting duplicates ###\n")
#  cat("###########################\n\n")

  pb <- txtProgressBar(min = 0, max = numfiles, style = 3)

  ucounts <- unlist(mclapply(1:numfiles, function(x){

    res<-(as.numeric(system(paste("sambamba markdup -r",bamFiles[x],"/dev/stdout 2> /dev/null | samtools view -c -"),intern=TRUE)))

    setTxtProgressBar(pb, x)

    return(res)

  }, mc.cores=cores, mc.preschedule=FALSE))

  close(pb)

  cat("### getting alignment stats ###\n")

  ascores=data.matrix(as.data.frame(
    mclapply(1:numfiles, function(x){


      res<-(as.numeric(readLines(pipe(paste(
        "samtools view",bamFiles[x],"| awk '{if($5==0){unaligned+=1}; if($5>0 && $5<",uniquescore,"){multimapped+=1}; if($5>=",uniquescore,"){unique+=1}}; END{print unaligned; print multimapped; print unique; print NR}'"
      )))))

      setTxtProgressBar(pb, x)

      return(res)

    }, mc.cores=11, mc.preschedule=FALSE )
  ))

  close(pb)
  cat("\n")

  barplot(data.matrix(ascores[3:1,])/1000000,beside=F,las=3,cex.names=0.7,ylim=c(0,40),col=c("blue","orange","red"))

  pb <- txtProgressBar(min = 0, max = numfiles, style = 3)

  cscores=t(data.matrix(do.call(cbind,mclapply(1:numfiles, function(x){


      res<-read.delim(pipe(paste(
        "samtools view -q",uniquescore,bamFiles[x],"| awk '{ a[$3]+=1 } END {for (i in a) print i,a[i] }' OFS='\t'"
      )), row.names=1,stringsAsFactors=FALSE,header=FALSE)

      setTxtProgressBar(pb, x)

      return(res)

  }, mc.cores=11, mc.preschedule=FALSE ))))

  colnames(ascores)<-sampleNames

  close(pb)
  cat("\n")

  cpscores <- sweep(cscores,1,rowSums(cscores),"/")

  cpscores<-cpscores[,mixedorder(colnames(cpscores))]
  barplot(cpscores,beside=T,las=3)


  rownames(cscores)=sampleNames


  barplot(cpscores,beside=T,las=3)
  legend("topright",legend=rownames(cpscores))

  ascores[is.na(ascores)]<-0
  rownames(ascores)<-c("unmapped","multimapped","unique","total")
  colnames(ascores)<-basename(bamFiles)
  ascores<-t(ascores)
  ascores<-as.data.frame(ascores)
  ascores$duplicate<-ascores$total-ucounts


  if( !is.null(probes) & expandprobes>0 ){
    probefile=bedSlop(probes,genomefile,expandprobes,expandprobes)
  } else{
    probefile=probes
  }

  if( !is.null(probes) ){
    olwp <- unlist(mclapply(1:numfiles, function(x){
      setTxtProgressBar(pb, x)
      as.numeric(system(paste("samtools view -c -L",probefile,"-q",uniquescore,bamFiles[x]),intern=TRUE))
    }, mc.cores=cores, mc.preschedule=F))
    ascores$overlappingprobes <- olwp
  }




  pscores<-sweep(ascores,1,ascores$total,"/")


  #uorder <- order(props[1,],decreasing=TRUE)
  #props=props[,uorder]
  #scores=scores[,uorder]

  if(plot){
    layout(matrix(c(1,2),nrow=2),heights=c(1,2))
    par(mar=c(1,4,6,6),xpd=FALSE)
    barplot(data.matrix(t(ascores[,1:4])),axisnames=F,xlab=F,ylab="number of reads",las=3,cex.names=0.75,col=c("blue","orange","red"),main="read alignment  quality statistics")
    abline(h=0,lwd=2)
    par(xpd=TRUE)
    legend(1.2*length(bamFiles),max(unlist(ascores)),yjust=0,legend=c("unique","multimapped","unmapped"),fill=c("blue","orange","red"))

    par(mar=c(20,4,1,6),xpd=FALSE)
    barplot(t(pscores[,1:3]),names=basename(removeext(bamFiles)),ylab="proportion of reads",las=3,cex.names=0.75,col=c("blue","orange","red"))
  }

  out<-list(ascores,pscores)
  names(out)=c("ascores","pscores")
  return(out)



}
