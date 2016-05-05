bgRead <- function ( bgfiles , makematrix=TRUE , bgnames=NULL , enforceEquality=FALSE , scores=TRUE , threads=getOption("threads",1L) ){

	numbgs<-length(bgfiles)



	if(scores){

		cmdString <- paste("cut -f 4",bgfiles)

		bgl <- cmdRun(cmdString,intern=T,threads=threads)

		bgl <- lapply(bgl,as.numeric)

		if(length(unique(unlist(lapply(bgl,length)))) == 1){equallength=TRUE} else{equallength=FALSE}

		if(!equallength & makematrix & enforceEquality){stop("number of scores among files are not identical")}

		if(is.null(bgnames)){ names(bgl) <- basename(removeext(bgfiles)) } else{ names(bgl) <- bgnames }

		if(makematrix & equallength){

			bgl<-data.matrix(as.data.frame(bgl))
			
		}

	} else{

		bgl <- mclapply( bgfiles , read_tsv , col_names=FALSE , mc.cores=threads, mc.preschedule=FALSE )
		if(is.null(bgnames)){ names(bgl) <- basename(removeext(bgfiles)) } else{ names(bgl) <- bgnames }

	}


	return(bgl)

}
