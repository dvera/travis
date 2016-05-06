bgQuantileNorm <-
function( bgfiles, normalizeto = 1:length(bgfiles) , threads=getOption("threads",1L) ){

	# check of files have identical coordinates in head and # of total lines
	# check if more than 1 file in bgfiles
	# make sure all files only have 4 columns

	fl=filelines(bgfiles,threads=threads)
	if(length(unique(fl))!=1){cat("WARNING: FILES HAVE UNEQUAL LINE COUNTS\n"); unequal=TRUE}

	numfiles<-length(bgfiles)
	bgnames<-basename(removeext(bgfiles))
	outnames<-paste0(bgnames,"_qnorm.bg")

	#cat("sorting reference\n")
	if(!unequal){
		bgrefs <- bgSort(bgfiles[normalizeto], threads=threads)
	} else{
		bgrefs=bgfiles[normalizeto]
	}

	if(length(normalizeto)>1){
		#stopifnot(length(unique(fl[normalizeto]))!=1)
		if(unequal){
			bgref <- bedCat(bgrefs)
		} else{
			bgref <- bgOps(bgrefs, operation="mean", outnames = paste0("mean_",basename(bgrefs[1])))
		}
	} else{
		bgref <- bgrefs
	}

	bgrefl<-filelines(bgref)

	if(any(fl>bgrefl)){
		bgref <- bedCat(rep(bgref,ceiling(fl/bgrefl)))
	}



	# if(unequal){
		cmdString <- paste("bash -c 'paste <(sort -T . -k4,4n",bgfiles,") <(shuf -n",fl,bgref," | sort -T . -k4,4n)' | awk '{print $1,$2,$3,$8}' OFS='\t' | sort -T . -k1,1 -k2,2n >",outnames)
		#cmdString <- paste0("paste <(sort -T . -k4,4n ",bgfiles,") <(shuf -n",fl,bgref," | sort -T . -k4,4n) | awk '{print $1,$2,$3,$8}' OFS='\t' | sort -T . -k1,1 -k2,2n >",outnames)
	# } else{
		# cmdString <- paste("sort -T . -k4,4n",bgfiles,"| paste -",bgref,"| awk '{print $1,$2,$3,$8}' OFS='\t' | sort -T . -k1,1 -k2,2n >",outnames)

	# }

	res <- cmdRun(cmdString,threads=threads)
	# unlink(bgref)
	# unlink(bgrefs)
	return(outnames)

}
