bgPoissonTest <- function( bedGraph1, bedGraph2=NULL, lower=FALSE, adjustp="fdr", threshold=0.05, plot=FALSE, lambda=NULL, trimmedmean=c(0.01,0.99), addone=TRUE, log10p=TRUE ){

  outname <- paste0(basename(removeext(bedGraph1)), "_ppois.bed")
  bgoutname <- paste0(basename(removeext(bedGraph1)), "_ppois.bg")

  if(is.null(bedGraph2)){
    paired=FALSE
  } else{
    paired=TRUE
    if(!is.null(lambda)){cat("bedGraph2 defined, ignoring lambda\n")}
  }

  bg1 <- read_tsv( bedGraph1, col_names=F)
  if(addone){ bg1[,4] <- bg1[,4]+1 }

  if(paired){
    bg2 <- read_tsv( bedGraph2, col_names=F)
    stopifnot(nrow(bg1)==nrow(bg2))
    if(addone){ bg2[,4] <- bg2[,4]+1 }
    lambdas <- bg2[,4]
    #lambdas[lambdas==0] <- NA
  } else if(!is.null(lambda)){
    lambdas <- lambda
    cat("using lambda=",lambdas,"\n")
  } else{
    lambdas <- mean(sort(bg1[,4])[round(trimmedmean[1]*length(bg1[,4])):round(trimmedmean[2]*length(bg1[,4])) ])
    cat("using lambda=",lambdas,"\n")
  }

  pvals <- ppois(bg1[,4],lambdas,lower=lower)
  #pvals[pvals==0]<-1

  if( is.null(adjustp) | is.na(adjustp) ){
    padjust=TRUE
    pvals2 <- pvals
  } else{
    pvals2 <- p.adjust(pvals, method=adjustp)
  }

  if( plot ){
    plot(density(pvals))
    lines(density(pvals2),col="red")
  }
  segs <- which(pvals2<=threshold)

  if(log10p){
    pvals2 <- -log10(pvals2)
    pvals2[is.infinite(pvals2)] <- 10
    pvals2[pvals2>10] <- 10
  }

  bgout <- bg1
  bgout[,4] <- pvals2
  write_tsv(bgout,bgoutname,col_names=FALSE)
  out <- cbind(bg1[segs,1:3],stringsAsFactors=F)
  write.tsv(out,file=outname)
  return(outname)
  #gaps <- c(bg1[segs[1],2],bg1[segs,3]) - c(bg1[segs,2],bg1[segs[length(segs)],3])




}
