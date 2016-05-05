bg.qnorm2 <-
function( bgs, numsamples=10, normalizeto=NULL, cores ="max" ){

	library(tools)
	library(parallel)

	numbgs<-length(bgs)

	if(cores=="max"){cores=detectCores()-1}
	if(cores > numbgs){cores=numbgs}


	bgnames<-basename(removeext(bgs))
	exts<-file_ext(bgs)
	outnames<-paste0(bgnames,"_qnorm.",exts)

	cat("reading\n")
	bglist<-mclapply(bgs,read.tsv,mc.cores=cores)
	numscores<-unlist(lapply(bglist,nrow))
	#bg<-read.tsv(bgs[1])

	cat("combining\n")
	if(is.null(normalizeto)){
		scores <- as.vector(unlist(lapply(bglist,"[",4)))
	} else{
		scores <- rep(as.vector(unlist(bglist[[normalizeto]][,4])),10)
	}

	cat("sampling\n")
	samples <- lapply(1:numbgs, function(y) sort(rowMeans(as.data.frame(mclapply(1:numsamples, function(x) sort(sample(scores,numscores[y])),mc.cores=cores)))) )

	cat("normalizing\n")
	bglist<-mclapply( 1:numbgs,function(x) {
			bglist[[x]][order(bglist[[x]][,4]),4]<-samples[[x]]
			bglist[[x]]
	}, mc.cores=cores )

	cat("saving\n")
	mclapply(1:numbgs, function(x) write.tsv(bglist[[x]],file=outnames[x]) ,mc.cores=cores )
	return(outnames)
}
