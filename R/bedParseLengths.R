#' Parse intervals in bed files by size.
#'
#' \code{bedParseLengths} separates intervals from a bed file into separate bed files based on interval sizes using the system \code{awk}.
#'
#' @param bedFiles A character vector of paths to bed files.
#' @param breaks A numeric vector of read lengths used as break points to parse reads by length. For example, if breaks=c(0,100,200), reads will be parsed into reads of length 0-99 and 100-199. Length of vector must beat least 2.
#' @param threads A positive integer specifying how many bams to process simultaneously.

bedParseLengths <-
function( bedFiles , breaks , threads=getOption("threads",1L) ){

	# check parameters
	numbreaks<-length(breaks)
	if(numbreaks==1){stop("need two or more breakpoints")}
	numfiles <- length(bedFiles)
	options(scipen=999999)


	# define output names
	bednames <- basename(removeext(bedFiles))
	brknames <- formatC(breaks,width=nchar(breaks[numbreaks]),format="d",flag="0")

	leftbreaks  <- brknames[ 1:(numbreaks-1) ]
	rightbreaks <- brknames[ 2:numbreaks     ]

	# parse fragments by length
	outnames<-lapply(1:numfiles,function(i){
		paste0(bednames[i],"_",leftbreaks,"-",rightbreaks,".bed")
	})

	cmdString<-unlist(lapply(1:numfiles, function(i){ paste (
			"awk '{ l=$3-$2; ",
			paste0("if(l >=",as.numeric(leftbreaks)," && l <",as.numeric(rightbreaks),"){ print $0 > \"",outnames[[i]],"\" }", collapse="" ),
			"}' OFS='\t'", bedFiles[i]
	)}))

	res <- cmdRun(cmdString,threads)

	return(outnames)
}
