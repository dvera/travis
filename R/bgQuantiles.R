bgQuantiles <-
function( bgs, quantiles=c(0,0.25,0.5,0.75,1), threads=getOption("threads",1L) ){


  numbgs<-length(bgs)



  bgnames<-basename(bgs)
  # exts<-file_ext(bgs)
  # outnames<-paste0(bgnames,"_iqrnorm.",exts)

  # count scores
  bglines <- filelines(bgs,threads=threads)
  ql <- data.matrix(as.data.frame(lapply(1:numbgs,function(x) round(bglines[x]*quantiles))))
  ql[ql==0] <- 1

  cmdString <- unlist(lapply(1:numbgs, function(x) {
    paste(
      "sort -T . -k4,4n", bgs[x],
      "| awk '{",
      paste("if (NR==",ql[,x],"){print $4}", collapse=";"),
      "}'"
    )
  }))

  q <- cmdRun( cmdString , threads=threads , lines=TRUE)
  q <- lapply(q,as.numeric)
  q <- data.matrix(as.data.frame(q))
  colnames(q) <- bgnames
  rownames(q) <- paste0("q",quantiles)

  return(q)
}
