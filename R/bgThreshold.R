bgThreshold <-
function( bgFiles, cutoff, mergeDist=50, positive=TRUE, negative=FALSE){

	bgnames<-basename(removeext(bgFiles))

	if(positive==TRUE){

		outnames<-paste0(bgnames,"_gt",abs(cutoff),".bg")

		cmdString <- paste(
			"awk '$4 >=",abs(cutoff),"'",bgFiles,
			if( !is.null(mergeDist) ){ "| bedtools merge -c 4 -o max -i stdin -d" },
			if( !is.null(mergeDist) ){ mergeDist },
			">",outnames
		)

		for(i in 1:length(bgFiles)){
			print( cmdString[i] )
			system( cmdString[i] )
		}
	}

	if(negative==TRUE){

		outnames<-paste0(bgnames,"_lt",-abs(cutoff),".bg")

		cmdString <- paste(
			"awk '$4 <=",-abs(cutoff),"'",bgFiles,
			if( !is.null(mergeDist) ){ "| bedtools merge -c 4 -o min -i stdin -d" },
			if( !is.null(mergeDist) ){ mergeDist },
			">",outnames
		)

		for(i in 1:length(bgFiles)){
			print( cmdString[i] )
			system( cmdString[i] )
		}
	}

	return(outnames)
}
