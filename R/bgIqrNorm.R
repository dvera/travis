bgIqrNorm <-
function( bgs, referenceIqr=NULL , mediancenter=FALSE, threads=getOption("threads",1L) ){


  numbgs<-length(bgs)
  bgnames<-basename(removeext(bgs))

  exts<-file_ext(bgs)
  outnames<-paste0(bgnames,"_iqrnorm.",exts)

  # calculate IQRs
  quantiles <- bgQuantiles ( bgs , threads=threads )
  iqrs <- quantiles[4,]-quantiles[2,]
  medians <- quantiles[3,]
  # use mean iqr as target IQR if referenceIqr not defined
  if(!is.null(referenceIqr) & is.numeric(referenceIqr) ){ iqr <- referenceIqr } else{ iqr <- mean(iqrs)}

  # calculate normalization factor
  scalars <- iqr/iqrs

  cmdString <- paste(
        "awk '{",
        if(mediancenter){"$4=$4-"},
        if(mediancenter){medians},
        if(mediancenter){";"},
        "$4=$4*",scalars,
        "; print $0}' OFS='\\t'",
        bgs,">",outnames
    )

  res <- cmdRun(cmdString, threads=threads)

  return(outnames)

}
