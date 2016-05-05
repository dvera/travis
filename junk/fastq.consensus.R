fastq.consensus <- function(fastqFiles, numSamples=1000, threshold=0.9, cores="max"){

  numfiles=length(fastqFiles)
  library(parallel)
  if(cores=="max"){cores=detectCores()-1}
  numSamples=rep(numSamples,numfiles)



  flist <- unlist(mclapply(1:numfiles,function(x){
    f=read.delim(pipe(paste("tail -n",numSamples[x]*4,fastqFiles[x],"| paste - - - - | cut -f 2")),stringsAsFactors=F,header=F)
    numSamples[x] <- nrow(f)
    fl<-t(as.data.frame(lapply(f,strsplit,split="")))
    colnames(fl)<-paste0("B",1:ncol(fl))
    rownames(fl)<-paste0("R",1:nrow(fl))
    fl=as.data.frame(fl,stringsAsFactors=F)
    fl[numSamples+1,]="A"
    fl[numSamples+2,]="T"
    fl[numSamples+3,]="C"
    fl[numSamples+4,]="G"
    fl[numSamples+5,]="N"
    flt1=as.data.frame(lapply(lapply(fl,table),as.vector),stringsAsFactors=F)-1
    rownames(flt1)<-c("A","C","G","N","T")
    flt2=as.data.frame((flt1-1)/numSamples > threshold)
    flt3=lapply(flt2,which)
    flt3[which(unlist(lapply(flt3,length))==0)]<-"N"
    flt3[which(unlist(lapply(flt3,unlist))==1)]<-"A"
    flt3[which(unlist(lapply(flt3,unlist))==2)]<-"C"
    flt3[which(unlist(lapply(flt3,unlist))==3)]<-"G"
    flt3[which(unlist(lapply(flt3,unlist))==4)]<-"N"
    flt3[which(unlist(lapply(flt3,unlist))==5)]<-"T"
    flt4<-paste(unlist(flt3),collapse="")
    return(flt4)
  }, mc.cores=cores))

  return(flist)

}
