bedParseStrands <-
function( bedfile , plus.suffix="plus" , minus.suffix="minus" ){

	library(tools)
	
	# make sure input is only 1 file
	if(length(bedfile) > 1){stop("bedParseLengths can only take 1 file")}

	ext<-file_ext(bedfile)
	fragname<-basename(removeext(bedfile))
	
	outplus <- paste0(fragname,"_",plus.suffix,".",ext)
	outminus <- paste0(fragname,"_",minus.suffix,".",ext)

	system(paste("awk '$6=\"+\"' OFS='\t'",bedfile,">",outplus))
	system(paste("awk '$6=\"-\"' OFS='\t'",bedfile,">",outminus))
	
	return(c(outplus,outminus))
}
