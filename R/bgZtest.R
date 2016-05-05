bgZtest <- function( testbg , controlbg , verbose=T , testname=NULL , controlname=NULL , windowsize=100 , pvalue=0.05 , padjmethod="fdr", threads=11 , mergewithin=50 ){

  # assumes dimensions of testbg and controlbg are the same and have identical coordinates

  if(is.null(testname)){ testname <- basename(removeext(testbg))}
  if(is.null(controlname)){ controlname <- basename(removeext(controlbg))}

  if(verbose){cat("reading in data\n")}
  bgt <- read_tsv(testbg,col_names=FALSE)
  bgc <- read_tsv(controlbg, col_names=FALSE)

  if(verbose){cat("finding unique windows\n")}
  #bgt2 <- aggregate(bgt[,4,drop=F],by=list(bgt[,1],bgt[,2],bgt[,3]),FUN="mean")
  bgt2 <- unique(bgt[,1:3])
  bgt2 <- bgt2[order(bgt2[,1],bgt2[,2],bgt2[,3]),]
  #bgc2 <- aggregate(bgc[,4,drop=F],by=list(bgc[,1],bgc[,2],bgc[,3]),FUN="mean")
  bgc2 <- unique(bgc[,1:3])
  bgc2 <- bgc2[order(bgc2[,1],bgc2[,2],bgc2[,3]),]

  numregions <- nrow(bgt)

  if(verbose){cat("splitting data by chromosome\n")}
  bgts <- split(bgt,bgt[,1])
  bgcs <- split(bgc,bgc[,1])
  bgts2 <- split(bgt2,bgt2[,1])
  bgcs2 <- split(bgc2,bgc2[,1])

  stopifnot(identical(unique(bgt[,1]),names(bgts)))

  bgtl <- unlist(lapply(bgts,nrow))

  if(verbose){cat("finding score windows\n")}
  nearby <- mclapply(seq_len(length(bgts)), function(r){ lapply(seq_len(nrow(bgts2[[r]])), function(i) { which( abs(bgts[[r]][,2] - bgts2[[r]][i,2]) <= windowsize ) } ) }, mc.cores=threads, mc.preschedule=FALSE )
  numnearby <- lapply(nearby,function(x) unlist(lapply(x,length)) )

  if(verbose){cat("calculating means\n")}
  tmeans <- lapply(seq_len(length(bgts2)), function(r){ unlist( lapply(seq_len(nrow(bgts2[[r]])), function(i) { mean(bgts[[r]][nearby[[r]][[i]],4] ) } ) ) } )
  cmeans <- lapply(seq_len(length(bgcs)), function(r){ unlist( lapply(seq_len(nrow(bgcs2[[r]])), function(i) { mean(bgcs[[r]][nearby[[r]][[i]],4] ) } ) ) } )

  #tsds <- lapply(seq_len(length(bgts)), function(r){ unlist( lapply(seq_len(nrow(bgts2[[r]])), function(i) { sd(bgts[[r]][nearby[[r]][[i]],4] ) } ) ) } )
  #csds <- lapply(seq_len(length(bgcs)), function(r){ unlist( lapply(seq_len(nrow(bgcs2[[r]])), function(i) { sd(bgcs[[r]][nearby[[r]][[i]],4] ) } ) ) } )
  if(verbose){cat("calculating variance\n")}
  tsd <- sd(bgt[,4])^2
  csd <- sd(bgc[,4])^2

  if(verbose){cat("calculating directionality of differences\n")}
  diffs <- unlist(lapply(1:length(bgts2), function(r){
    ( tmeans[[r]]-cmeans[[r]] )
  }))

  if(verbose){cat("calculating z-scores\n")}
  zvals <- unlist(lapply(1:length(bgts2), function(r){
    ( tmeans[[r]]-cmeans[[r]] ) / sqrt( tsd/numnearby[[r]] +  csd/numnearby[[r]] )
  }))

  if(verbose){cat("transforming p-values\n")}
  pvals <- pnorm( -abs(zvals),0,1 )

  bgp <- bgt2
  bgp[,4] <- pvals
  bgp[,5] <- diffs

  if(is.null(padjmethod) | is.na(padjmethod)) {
    adjustname <- "noAdj"
    bgp[,4] <- bgp[,4]
  } else{
    adjustname <- padjmethod
    bgp[,4] <- p.adjust( bgp[,4], padjmethod )
  }

  signs <- sign(diffs)
  signs[signs==0] <- 1
  bgp[,4] <- signs * -log10(bgp[,4])




  if(verbose){cat("thresholding\n")}
  bgseg <- bgp[which(abs(bgp[,4]) >= -log10(pvalue) ) ,]


  bgpb <- data.frame(bgseg[,1],bgseg[,2],bgseg[,3],paste0("seg",1:nrow(bgseg)),round(abs(100*bgseg[,4])),"+",stringsAsFactors=FALSE)
  bgpb[which(bgseg[,5]<0),6] <- "-"

  posseg <- bgpb[which(bgpb[,6]=="+"),1:3]
  negseg <- bgpb[which(bgpb[,6]=="-"),1:3]


  bgoutname <- paste0(testname,"_vs_",controlname,"_zTest",adjustname,".bg")
  possegoutname <-  paste0(testname,"_vs_",controlname,"_zTest",adjustname,"_pos_",as.character(100*pvalue),"pct.bed")
  negsegoutname <-  paste0(testname,"_vs_",controlname,"_zTest",adjustname,"_neg_",as.character(100*pvalue),"pct.bed")


  if(verbose){cat("saving data\n")}
  write.tsv(bgp[,1:4],file=bgoutname)
  write.tsv(posseg,file=possegoutname)
  write.tsv(negseg,file=negsegoutname)

  if(verbose){cat("conforming bg\n")}
  bgsegc <- bgConform(bgoutname,threads=threads)

  if(verbose){cat("merging segments\n")}
  possegm <- bedtoolsMerge(possegoutname,flank=mergewithin)
  negsegm <- bedtoolsMerge(negsegoutname,flank=mergewithin)

  return(list(bg=bgsegc,beds=c(possegm,negsegm)))

}

bgConform <- function( bgFiles, threads=getOption("threads",1L) ){

  outnames <- paste0(basename(removeext(bgFiles)),"_conform.bg")

  cmdString <- paste(
    #"cut -f 1,2,3",bgFiles,
    #"| bedtools merge -i /dev/stdin",
    "bedtools merge -i", bgFiles,
    "| bedtools makewindows -b /dev/stdin -w 1",
    "| bedtools map -b", bgFiles,"-a /dev/stdin",
    "-o mean -c 4",
    #if(reduce){"| awk '{if(NR==1){a=$1;b=$2;c=$3;d=$4; }'"},
    ">", outnames
  )

  cmdRun( cmdString, threads=threads)

  return(outnames)

}
