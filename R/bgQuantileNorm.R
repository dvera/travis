bgQuantileNorm <-
function( bgfiles, normalizeto = 1:length(bgfiles) , threads=getOption("threads",1L) ){

	# check of files have identical coordinates in head and # of total lines
	# check if more than 1 file in bgfiles
	# make sure all files only have 4 columns

	numfiles<-length(bgfiles)
	bgnames<-basename(removeext(bgfiles))
	outnames<-paste0(bgnames,"_qnorm.bg")

	#cat("sorting reference\n")
	bgrefs <- bgSort(bgfiles[normalizeto], threads=threads)

	if(length(normalizeto)>1){
		bgref <- bgOps(bgrefs, operation="mean", outnames = paste0("mean_",basename(bgrefs[1])))
	} else{
		bgref <- bgrefs
	}

	cat("normalizing data\n")
	cmdString <- paste("sort -T . -k4,4n",bgfiles,"| paste -",bgref,"| awk '{print $1,$2,$3,$8}' OFS='\t' | sort -T . -k1,1 -k2,2n >",outnames)

	res <- cmdRun(cmdString,threads=threads)
	unlink(bgref)
	unlink(bgrefs)
	return(outnames)

}
