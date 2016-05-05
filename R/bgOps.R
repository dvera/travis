bgOps <-
function( bglist1 , operation , bglist2=NULL , outnames = NULL , pattern=NULL , replacement=NULL , threads=getOption("threads",1L), forceall=FALSE ){

	library(tools)

	# check that file has only 4 columns
	# check that file line counts are equal

	numbgs <- length(bglist1)

	############################################
	############# argument check ###############
	############################################
	if( is.null(outnames)==FALSE & is.null(pattern)==FALSE & is.null(replacement)==FALSE & numbgs!=1 & is.null(bglist2)==FALSE) {
		stop("must use pattern/replacement OR outnames, not both")
	}

	if( length ( operation ) != 1){stop("can only perform 1 operation")}

	if(is.null(bglist2)==FALSE & length(bglist2) == 1 & numbgs > 1){
		bglist2 <- rep(bglist2,numbgs)
	}

	bgnames1=basename(removeext(bglist1))
	if(!is.null(bglist2)){bgnames2=basename(removeext(bglist2))}

	if( is.null(outnames) ){
		outnames<-gsub(pattern,replacement,basename(bglist1))

		if(length(which(bglist1 %in% outnames)==TRUE)>0){stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}

		if(is.null(bglist2)==FALSE){
			if(length(which(bglist2 %in% outnames)==TRUE)>0) {stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}
		}
	}

	if( is.null(bglist2)==FALSE & length(outnames) != numbgs){ stop("outnames must be of same length as bglist1")}

	############################################
	############################################
	############################################

	if( is.null(bglist2) & operation %in% c("mean","sd")){

		# define pasted score columns
		scorecols<-paste(paste0("$",4*(1:numbgs)),collapse="+")

		# make sure outnames is length 1
		if(is.null(outnames)){stop("must define a single file name in 'outnames' for mean/sd operations on a set of files")}

		# perform operation
		if(operation == "mean" ){
			cmdString <- paste0("paste ",paste(bglist1,collapse=" ")," | awk '{print $1,$2,$3,(",scorecols,")/",numbgs,"}' OFS='\t' > ",outnames)
		} else if(operation == "sd"){
			cmdString <- paste0("paste ",paste(bglist1,collapse=" ")," | awk 'BEGIN {n=",numbgs,"} { s=0 ; m=(",scorecols,")/",numbgs,"; for (i=1;i<=",numbgs,";i++){ s+=($(i*4)-m)^2 } print $1,$2,$3,sqrt(s/(n-1))}' OFS='\t' > ",outnames)
		}

	}

	if( is.null(bglist2) & grepl("\\$",operation) ){
		cmdString <- paste0("awk '{print $1,$2,$3,",operation,"}' OFS='\t' ",bglist1," > ",outnames)
	}


	if( is.null(bglist2) & operation %in% c("log2","log10","antilog2","antilog10","inverse","mediancenter","meancenter")){

		# check operation
		if(operation %ni% c("mean","sd","log2","log10","antilog2","antilog10","inverse","mediancenter","meancenter")) { stop("operations that can be performed on a single set of bedGraphs include mean, log2, log10, antilog2, antilog10, mediancenter, meancenter, and inverse")}

		# define output file names
		if(is.null(outnames)){outnames <- paste0 (bgnames1 , "_",operation,".bg")}

		if(operation=="log2"){
			cmdString <- paste("awk '{print $1,$2,$3,log($4)/log(2)}' OFS='\t' ",bglist1," > ",outnames)
		}
		if(operation=="log10"){
			cmdString <- paste("awk '{print $1,$2,$3,log($4)/log(10)}' OFS='\t' ",bglist1," > ",outnames)
		}
		if(operation=="antilog2"){
			cmdString <- paste("awk '{print $1,$2,$3,2^$4}' OFS='\t' ",bglist1," > ",outnames)
		}
		if(operation=="antilog10"){
			cmdString <- paste("awk '{print $1,$2,$3,10^$4}' OFS='\t' ",bglist1," > ",outnames)
		}
		if(operation=="inverse"){
			cmdString <- paste("awk '{print $1,$2,$3,1/$4}' OFS='\t' ",bglist1," > ",outnames)
		}
		if(operation=="mediancenter"){
			cmdString <- paste("SCOREMEAN=$(cut -f 4",bglist1,"| sort -k1,1n | awk 'BEGIN{i=0} {a[i]=$1 ; i++} END{print a[int(NR/2)]}') && awk -v scoremean=\"$SCOREMEAN\" '{print $1,$2,$3,$4-scoremean}' OFS='\t'",bglist1,">",outnames)
		}
		if(operation=="meancenter"){
			cmdString <- paste("SCOREMEAN=$(awk '{a+=$4} END {print a/NR}'",bglist1,") && awk -v scoremean=\"$SCOREMEAN\" '{print $1,$2,$3,$4-scoremean}' OFS='\t' ",bglist1," > ",outnames)
		}


	}

	if( is.null(bglist2)==FALSE & length(bglist2) == numbgs ){

		# check operation
		if(operation %ni% c("log2ratio","ratio","difference","mean","absdiff","sum")) { stop("operations that can be performed on pairs of files include log2ratio, ratio, difference, and mean")}

		# define output file names
		if(is.null(outnames)){outnames <- paste0 (bgnames1 , "_",operation,"_",bgnames2,".bg")}

		# perform operation
		if(operation=="log2ratio"){
			cmdString <- paste("paste",bglist1,bglist2,"| awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}",if(forceall){"else{print $1,$2,$3,0}"},"}' OFS='\t' >",outnames)
		}
		if(operation=="ratio"){
			cmdString <- paste("paste",bglist1,bglist2,"| awk '{if($8 != 0){print $1,$2,$3,($4/$8)}",if(forceall){"else{print $1,$2,$3,0}"},"}' OFS='\t' >",outnames)
			#cmdString <- paste("paste",bglist1,bglist2,"| awk '{print $1,$2,$3,$4/$8}' OFS='\t' >",outnames)
		}
		if(operation=="difference"){
			cmdString <- paste("paste",bglist1,bglist2,"| awk '{print $1,$2,$3,$4-$8}' OFS='\t' >",outnames)
		}
		if(operation=="absdiff"){
			cmdString <- paste("paste",bglist1,bglist2,"| awk '{a=$4-$8; b=sqrt(a*a); print $1,$2,$3,b}' OFS='\t' >",outnames)
		}
		if(operation=="mean"){
			cmdString <- paste("paste",bglist1,bglist2,"| awk '{print $1,$2,$3,($4+$8)/2}' OFS='\t' >",outnames)
		}
		if(operation=="sum"){
			cmdString <- paste("paste",bglist1,bglist2,"| awk '{print $1,$2,$3,$4+$8}' OFS='\t' >",outnames)
		}


	}

	res <- cmdRun( cmdString, threads )

	return(outnames)

}
