bedEnds <-
function( bedFiles, end=5, strand=TRUE , threads=getOption("threads",1L) ){

	outnames<-paste0(basename(removeext(bedFiles)),"_",if(end==5){"5p"}else{"3p"},if(strand){"s"},".bed" )

	cmdString <- paste(
		if(strand){
			if(end==5){
				"awk '{if ($6==\"-\"){ $2=$3-1 } else{ $3=$2+1 } ; if($7){$7=$2; $8=$3} ; print $0}'"
			}
			if(end==3){
				"awk '{if ($6==\"-\"){ $3=$2+1 } else{ $2=$3-1 } ; if($7){$7=$2; $8=$3} ; print $0}'"
			}
		} else{
			if(end==5){
				"awk '{ $2=$3-1 } ; if($7){$7=$2; $8=$3} ; print $0 }'"
			}
			if(end==3){
				"awk '{ $3=$2+1 } ; if($7){$7=$2; $8=$3} ; print $0 }'"
			}
		}, "OFS='\t'", bedFiles,
		"| sort -T . -k1,1 -k2,2n >",
		outnames )

	cmdRun(cmdString,threads)

	return(outnames)

}
