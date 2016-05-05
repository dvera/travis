bedCenters <-
function( bedfiles, postsort=TRUE , buffer="10G", sizeScore=FALSE, threads=getOption("threads",1L) ){

	#MAKE A FILE OF FRAGMENT CENTERS

	outnames<-paste(basename(removeext(bedfiles)),"_centers.",if(sizeScore){"bg"} else{"bed"},sep="")

	cmdString <- paste0(
		"awk '{a=int(($2+$3)/2+0.5);",if(sizeScore){"$4=$3-$2;"},"$2=a; $3=a+1;print}' OFS='\t' ",
		bedfiles,
		if(postsort){paste(" | sort -T . -S ",buffer," -k1,1 -k2,2n")},
		">",outnames
	)

	res <- cmdRun(cmdString,threads)

	return(outnames)

}
