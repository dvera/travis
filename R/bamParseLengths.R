#' Parse reads in bam files by insert size.
#'
#' \code{bam.parselengths} separates reads from a paired-end bam file into separate files based on insert sizes using \code{samtools view} and the system \code{awk}.
#'
#' @param bamFiles A character vector of paths to bam files.
#' @param breaks A numeric vector of read lengths used as break points to parse reads by length. For example, if breaks=c(0,100,200), reads will be parsed into reads of length 0-99 and 100-199. Length of vector must beat least 2.
#' @param threads A positive integer specifying how many bams to process simultaneously.

bamParseLengths <-
function( bamFiles , breaks , threads=getOption("threads",1L) ){

	# check parameters
	numbreaks<-length(breaks)
	if(numbreaks==1){stop("need two or more breakpoints")}
	numfiles <- length(bamFiles)
	options(scipen=999999)


	# define output names
	bamnames <- basename(removeext(bamFiles))
	brknames <- formatC(breaks,width=nchar(breaks[numbreaks]),format="d",flag="0")

	leftbreaks  <- brknames[ 1:(numbreaks-1) ]
	rightbreaks <- brknames[ 2:numbreaks     ]

	# parse fragments by length
	outnames<-lapply(1:numfiles,function(i){
		paste0(bamnames[i],"_",leftbreaks,"-",rightbreaks,".sam")
	})

	cmdString<-unlist(lapply(1:numfiles, function(i){
		paste (
			"samtools view -h",bamFiles[i],"| awk '{if($1~\"@\"){",
				paste0("print $0 > \"", outnames[[i]],"\"",collapse=";"),
			"}",
			paste0("else if(sqrt($9*$9) >= ",leftbreaks," && sqrt($9*$9) < ",rightbreaks,"){ print $0 > \"",outnames[[i]],"\" }", collapse="" ),
			"}' OFS='\t'"
		)
	}))

	res <- cmdRun(cmdString,threads)

	return(outnames)
}
