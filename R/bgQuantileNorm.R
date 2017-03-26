bgQuantileNorm <-
function( bgfiles, normalizeto = 1:length(bgfiles) , threads=getOption("threads",1L) ){

	# check of files have identical coordinates in head and # of total lines
	# check if more than 1 file in bgfiles
	# make sure all files only have 4 columns

	fl=filelines(bgfiles,threads=threads)
	if(length(unique(fl))!=1){cat("WARNING: FILES HAVE UNEQUAL LINE COUNTS\n"); unequal=TRUE} else{unequal=FALSE}

	numfiles<-length(bgfiles)
	bgnames<-basename(removeext(bgfiles))
	outnames<-paste0(bgnames,"_qnorm.bg")

	#cat("sorting reference\n")
	if(!unequal){
		cat("input data sets have equal lines\n")
		bgrefs <- bgSort(bgfiles[normalizeto], threads=threads)

	} else{
		cat("input data sets have unequal lines\n")
		bgrefs=bgfiles[normalizeto]
	}

	if(length(normalizeto)>1){
		#stopifnot(length(unique(fl[normalizeto]))!=1)
		if(unequal){
			lineratio <- ceiling(max(fl)/sum(fl[normalizeto]))
			tmpname <- randomStrings()
			tmpname2 <- randomStrings()
			bgref <- bedCat(bgrefs,tmpname)
			if(lineratio>1){
				bgref <- bedCat(rep(bgref,lineratio),tmpname2)
			}
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



	if(unequal){
		cmdString <- paste("bash -c 'paste <(sort -g -T . -k4,4n",bgfiles,") <(shuf -n",fl,bgref," | sort -g -T . -k4,4n)' | awk '{print $1,$2,$3,$8}' OFS='\t' | sort -T . -k1,1 -k2,2n >",outnames)
		#cmdString <- paste0("paste <(sort -T . -k4,4n ",bgfiles,") <(shuf -n",fl,bgref," | sort -T . -k4,4n) | awk '{print $1,$2,$3,$8}' OFS='\t' | sort -T . -k1,1 -k2,2n >",outnames)
	} else{
		bgref <- bgSort(bgref)
		cmdString <- paste("bash -c 'paste <(sort -g -T . -k4,4n",bgfiles,") <(sort -g -T . -k4,4n",bgref,")' | awk '{print $1,$2,$3,$8}' OFS='\t' | sort -T . -k1,1 -k2,3n >",outnames)

	}

	res <- cmdRun(cmdString,threads=threads)
	# unlink(bgref)
	# unlink(bgrefs)
	return(outnames)

}
