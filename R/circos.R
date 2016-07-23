circosPlot <- function( trackdata , types ,  windowsize=1000000 , widths=40 , colors="black" , cutoffs="n" , scales=TRUE, backgrounds=TRUE , threads=getOption("threads",1L) ){

  numtracks=length(trackdata)


  if(is.vector(colors)){colors=as.list(colors)}

  if(lengths(widths)==1){
    widths=rep(widths,numtracks)
  } else{
    if(length(widths)!=length(trackdata)){
      stop("lengths of widths should equal length of tracks")
    }
  }

  if(is.character(cutoffs) & lengths(cutoffs)==1){
    cutoffs=rep(cutoffs,numtracks)
  } else{
    if(length(cutoffs)!=length(trackdata)){
      stop("lengths of cutoffs should equal length of tracks")
    }
    if(is.vector(cutoffs)){
      cutoffs <- as.list(cutoffs)
    }
  }


  if(lengths(scales)==1){
    scales=rep(scales,numtracks)
  } else{
    if(length(scales)!=length(trackdata)){
      stop("lengths of scales should equal length of tracks")
    }
  }

  if(lengths(backgrounds)==1){
    backgrounds=rep(backgrounds,numtracks)
  } else{
    if(length(backgrounds)!=length(trackdata)){
      stop("lengths of backgrounds should equal length of tracks")
    }
  }

  datalist <- lapply(tsvRead(trackdata,threads=threads),as.data.frame)

  csi <- grep("chrom.sizes",trackdata)

  if(length(csi)==0){stop("can't find chrom sizes in trackdata")}
  chromsizes <- trackdata[csi[1]]

  chrombed <- as.data.frame(tsvRead(bedtoolsMakeWindows(chromsizes,windowsize,genome=T)))
  chrombed$X4=NA
  chrombed$X5=NA
  colnames(chrombed)=c("seg.name","seg.Start","seg.End","the.v","NO")


  db=segAnglePo(chrombed,seg=unique(chrombed[,1]))

  datalist[[csi[1]]] <- db

  rs=c(400,400-cumsum(widths))

  par(mar=c(1,1,1,1))
  plot(c(1,800),c(1,800),type='n',axes=F,xlab="",ylab="")

  hmps=c("topleft","topright","bottomright","bottomleft",rep("far",100))
  hmp=1

  for(i in seq_len(numtracks)){
    print(i)
    if(i==csi){
      circos(R=rs[i],cir=db,W=widths[i],type=types[i],scale=scales[i],B=backgrounds[i],col=colors[[i]],print.chr.lab = TRUE)
    } else{
      circos(R=rs[i],cir=db,W=widths[i],mapping=datalist[[i]],type=types[i],scale=scales[i],B=backgrounds[i],col=colors[[i]],col.v=4,col.bar=TRUE, col.bar.po=hmps[hmp],cutoff=cutoffs[[i]])
      if(grepl("heatmap",types[[i]])){hmp=hmp+1}
    }
  }



}
