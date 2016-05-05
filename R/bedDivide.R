bedDivide <- function( bed,regionsize,windowsize,covbedname, flank=1000, start=2,stop=3 , meta=FALSE ){

  options(scipen=99999)

  if(meta){

    curbed<-read_tsv(bed,col_names=FALSE)

    bedrows<-nrow(curbed)
    bedcols<-ncol(curbed)

    numwindows<-regionsize/windowsize
    numflankwindows<-flank/windowsize
    leftwinstarts<-0:(numflankwindows-1) * windowsize - flank
    leftwinends<-1:numflankwindows * windowsize - flank
    rightwinstarts<-0:(numflankwindows-1) * windowsize
    rightwinends<-1:numflankwindows * windowsize
    genesizes<-curbed[,stop]-curbed[,start]
    sizecos<-genesizes/regionsize
    genewinstarts<-mclapply(1:bedrows, function(x) curbed[x,start] + 0:(numwindows-1) * windowsize * sizecos[x], mc.cores=detectCores())
    genewinends<-mclapply(1:bedrows, function(x) curbed[x,start] + 1:numwindows * windowsize * sizecos[x], mc.cores=detectCores())
    genewinstarts<-mclapply(genewinstarts,round,mc.cores=detectCores())
    genewinends<-mclapply(genewinends,round,mc.cores=detectCores())

    covbed<-data.frame(
      "V1"=rep(curbed[,1],each=numwindows+numflankwindows*2),
      "V2"=format(unlist(lapply(1:bedrows,function(x){ c(leftwinstarts+curbed[x,start],genewinstarts[[x]],rightwinstarts+curbed[x,stop]) } )),scientific=F,trim=T),
      "V3"=format(unlist(lapply(1:bedrows,function(x){ c(leftwinends+curbed[x,start], genewinends[[x]],rightwinends+curbed[x,stop]) } )),scientific=F,trim=T),
      "V4"=1:(bedrows*(numwindows+numflankwindows*2))
    )

    write.tsv(covbed,file=covbedname)
    system(paste("sort -k1,1 -k2,2n -k3,3n",covbedname,"-o",covbedname))
    covbedorder.call <- pipe(paste("cut -f 4",covbedname),open="r")
    covbedorder<-as.numeric(readLines(covbedorder.call))
    close(covbedorder.call)
    return(covbedorder)

  } else{

    numwindows<-(regionsize/windowsize)

    bedname<-basename(removeext(bed))
    curbed<-read.tsv(bed)

    bedrows<-nrow(curbed)
    bedcols<-ncol(curbed)
    #check parameters
    if(ceiling(numwindows)!=floor(numwindows)){stop("regionsize is not a multiple of windowsize")}
    if(ceiling(regionsize)!=floor(regionsize)){stop("regionsize must be an even number")}

    #make covbed
    flanks<- (  0:(numwindows-1)  ) * windowsize
    winstarts<-as.numeric( unlist( lapply( curbed[,2] , function(x) x + flanks ) ) )
    covbed<-data.frame(
      "V1"=rep(curbed[,1],each=numwindows),
      "V2"=format(winstarts,scientific=F,trim=T),
      "V3"=format(winstarts+windowsize,scientific=F,trim=T),
      "V4"=1:(bedrows*numwindows)
    )
    write.tsv(covbed,file=covbedname)
    system(paste("sort -k1,1 -k2,2n -k3,3n",covbedname,"-o",covbedname))
    covbedorder.call <- pipe(paste("cut -f 4",covbedname),open="r")
    covbedorder<-as.numeric(readLines(covbedorder.call))
    close(covbedorder.call)
    return(covbedorder)
}
}
