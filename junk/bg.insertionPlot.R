bg.insertionPlot <- function( siteBgs , insertBgs , siteLoc="chr14:93749062" , insertLoc=c("chr16:48280951-48485357", "chr6:36534887-36740318","chr6:36534887-36740318") , sspan=0.1 , ispan=0.5 , siteColors=NULL , plotFlank=1000000 , cores="max" , ylims=c(-5,1) , insertNames=basename(removeext(insertBgs)) , siteNames=basename(removeext(siteBgs)) , linewidth=3 , i2s=NULL , ylabel="RT" , maintitle="" , linetypes=1 ){

  options(scipen=99999)
  library(parallel)
  if(cores=="max"){cores=detectCores()-1}

  #insertBgs <- siteBgs[i2s]

  # count files
  numSiteBgs <- length(siteBgs)
  numInsertBgs <- length(insertBgs)

  # if(length(sspan)==1){sspans=rep(sspan,numSiteBgs)} else{sspans<-sspan}
  # if(length(sspans) != numSiteBgs){stop("must have same number of sspans and siteBgs")}
  #
  # if(length(ispan)==1){ispan=rep(ispan,numInsertBgs)} else{ispans<-ispan}
  # if(length(ispan) != numInsertBgs){stop("must have same number of sspans and siteBgs")}

  if(length(linetypes)==1){linetype=rep(linetypes,numSiteBgs)} else{linetype<-linetypes}
  if(length(linetype) != numSiteBgs){stop("must have same number of sspans and siteBgs")}

  # set plot colors
  if(is.null(siteColors)){(siteColors=rainbow(numSiteBgs))}
  #if(is.null(insertColors)){(insertColors=rainbow(numInsertBgs))}
  #insertColors=siteColors[i2s]

  # extract site location
  siteLoc=gsub(",","",siteLoc)
  schrm <- unlist(strsplit(siteLoc,":"))[1]
  sstar <- as.numeric(unlist(strsplit(siteLoc,":"))[2])

  # extract insert location
  insertLoc=gsub(",","",insertLoc)
  ichrm <- unlist(lapply(strsplit(insertLoc,":"),"[",1))
  istar <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(insertLoc,":"),"[",2 ) ),"-"),"[",1)))
  istop <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(insertLoc,":"),"[",2 ) ),"-"),"[",2)))

  # extract insert data
  ibgl <- mclapply(1:numInsertBgs, function(x){
    read.delim(pipe(paste0("echo -e '",ichrm[x],"\t",istar[x],"\t",istop[x],"' | bedtools intersect -sorted -u -a ",insertBgs[x]," -b stdin")),stringsAsFactors=F,header=F)
  },mc.cores=cores)

  # instract site data
  sbgl <- mclapply(siteBgs, function(x){
    read.delim(pipe(paste0("echo -e '",schrm,"\t",sstar-(plotFlank*1.5),"\t",sstar+(plotFlank*1.5),"' | bedtools intersect -sorted -u -a ",x," -b stdin")),stringsAsFactors=F,header=F)
  },mc.cores=cores)


  # extract, smooth, and break insert data
  icoords <- lapply(lapply(ibgl,"[",2),unlist)
  icoords <- lapply(1:length(icoords),function(x){
    icoords[[x]]-(icoords[[x]][1])
  })
  icoords <- lapply(1:length(icoords),function(x){
    ic <- icoords[[x]]
    ic <- ic*(max(unlist(icoords)/max(ic)))
    ic <- ic+sstar+1
    return(ic)
  })




  # extract site data
  sscores <- as.data.frame(lapply(lapply(sbgl,"[",4),unlist))
  scoords <- sbgl[[1]][,2]
  sscores<-data.frame(lapply(1:numSiteBgs,function(x){

    if(sspan>0){
      loess(sscores[,x]~scoords,span=sspan)$fitted
    } else{
      sscores[,x]
    }

  }))
  # get insert size
  allicoords <- unlist(icoords,"[",2)
  isize <- max(allicoords)-min(allicoords)

  # extract site target
  slength<-length(scoords)
  sbreak <- which(scoords>sstar)[1]

  # break site target
  gap<-as.data.frame(matrix(NA,nrow=1,ncol=numSiteBgs))
  colnames(sscores)<-paste0("V",1:ncol(sscores))
  colnames(gap)<-paste0("V",1:ncol(sscores))
  sscores<-rbind(data.matrix(sscores[1:(sbreak-1),]),gap,data.matrix(sscores[sbreak:slength,]))

  scoords<-c(scoords[1:(sbreak-1)],scoords[sbreak-1]+1,scoords[sbreak:slength])
  scoords[(sbreak+1):(slength+1)]<-scoords[(sbreak+1):(slength+1)]+isize
  #sscores<-rbind(sscores[1:(sbreak-1),],gap,sscores[sbreak:slength,])

  par(mfrow=c(1,2))

  # plot data
  plot(0,type="n",ylim=ylims,xlim=c(sstar-plotFlank,sstar+isize+plotFlank),ylab=ylabel,xlab=paste0(schrm," coordinate (bp)"),main="smoothing before inserting")

  abline(v=c(sstar,sstar+isize),h=0,col="grey50")

  legend("topright",legend=siteNames,lwd=linewidth,col=siteColors,cex=0.7,lty=linetypes)
  #legend("top",legend=insertNames,lwd=linewidth,col=insertColors,cex=0.7)

  for(i in 1:numSiteBgs){
    lines(scoords,sscores[,i],col=siteColors[i],lwd=linewidth,lty=linetypes[i])
  }

  iscores <- lapply(lapply(ibgl,"[",4),unlist)
  iscores2<-iscores
  iscores<-lapply(1:numInsertBgs,function(x){
    if(ispan>0){loess(iscores[[x]]~icoords[[x]],span=ispan)$fitted  } else{iscores[[x]]}
  } )

  for(i in 1:numInsertBgs){
    lines(icoords[[i]],iscores[[i]],col=siteColors[i2s[i]],lwd=linewidth,lty=linetypes[i2s[i]])
  }



  if(!is.null(i2s)){
    plot(0,type="n",ylim=ylims,xlim=c(sstar-plotFlank,sstar+isize+plotFlank),ylab=ylabel,xlab=paste0(schrm," coordinate (bp)"),main="smoothing after inserting")

    abline(v=c(sstar,sstar+isize),h=0,col="grey50")

    legend("topright",legend=siteNames,lwd=linewidth,col=siteColors,cex=0.7,lty=linetypes)
    #legend("top",legend=insertNames,lwd=linewidth,col=insertColors)

    for(i in (1:numSiteBgs)[-i2s]){
      lines(scoords,sscores[,i],col=siteColors[i],lwd=linewidth,lty=linetypes[i])
    }


    for(i in 1:numInsertBgs){
      acoords <- data.frame(c(icoords[[i]],scoords),c(iscores2[[i]],sscores[,i2s[i]]))
      acoords <- acoords[order(acoords[,1]),]
      acoords <- acoords[complete.cases(acoords),]
      if(sspan>0){
        acoords[,2]<-loess(acoords[,2]~acoords[,1],span=sspan)$fitted
      } else{
        acoords[,2]
      }

      lines(acoords[,1],acoords[,2],col=siteColors[i2s[i]],lwd=linewidth,lty=linetypes[i])
    }

  }


}
