bed.recenter2 <-
function( bed,regionsize,center=FALSE,strand=TRUE,start=2,stop=3){
	
	bedname<-basename(removeext(bed))
	windowbed<-read.tsv(bed)
	outname<-paste0(bedname,"_recenter.bed")
	if( ceiliing(regionsize/2) != floor(regionsize/2) ){stop("regionsize must be divisible by 2")}
	# system(paste0("awk '{ ",if(strand=="TRUE"){"$2=int(($3-$2)/2)-",regionsize/2,
	# 	"if else($6==\"-\"){ $2=$",
	# 	stop,"-",regionsize/2,
	# 	"} else{ $2=$",start,"-",regionsize/2,
	# 	"}; $3=$2+",regionsize,
	# 	"; if($2>0){ print $0 } }' OFS='\t' ",
	# 	bed," | sort -k1,1 -k2,2n >",outname
	# ))

	return(outname)
}

