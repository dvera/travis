fastq.qc <- function(fastqFiles, threshold=20, numBadBases=10 , cores="max"){

  if(cores=="max"){cores=detectCores()-1}
  numfiles<-length(fastqFiles)
  if(cores>numfiles){cores<-numfiles}
  scorelist<-mclapply(1:numfiles, function(x){
    scores<-readLines(pipe(paste("cat",fastqFiles[x],"| paste - - - - | cut -f 4")))
    scores<-lapply(scores,strsplit,split="")
    scores<-lapply(scores,unlist)
    scores<-lapply(scores,lapply,charToRaw)
    scores<-lapply(scores,unlist)
    scores<-lapply(scores,as.numeric)
    scores<-lapply(1:length(scores), function(x){ scores[[x]]-33})
    return(scores)
  } , mc.cores=cores)



  lengths<-lapply(scorelist,lapply,length)
  lengths<-lapply(lengths,unlist)
  maxLength<-max(unlist(lapply(lengths,max)))

  rb<-rainbow(numfiles)
  plot(0,type="n",xlim=c(0,maxLength),ylim=c(0,1),xlab="cycle number",ylab="density of read lengths")
  for(i in 1:numfiles){
    lines(density(lengths[[i]]),col=rb[i])
  }


  cycleList<-mclapply(1:numfiles, function(x) lapply(1:maxLength,function(y) unlist(lapply(scorelist[[x]],"[",y))), mc.cores=cores)
  cycleList<-mclapply(cycleList,lapply,na.omit, mc.cores=cores)

  #cycleCoverage<-lapply(1:numfiles, function(x) lapply(1:maxLength, function(y) length(cycleList[[x]][[y]])))
  cycleCoverage<-unlist(lapply(cycleList,lapply,length))

  cyclePass<-mclapply(1:numfiles, function(x) unlist(lapply(1:maxLength, function(y) length(which(cycleList[[x]][[y]]>=20)))), mc.cores=cores)






  readPass<-unlist(mclapply(1:numfiles, function(x) length( which ( unlist(lapply(1:length(scorelist[[x]]), function(y) length(which(scorelist[[x]][[y]]<20 )))) > numBadBases ) ),mc.cores=cores))

  readCounts<-unlist(lapply(scorelist,length))
  
}
