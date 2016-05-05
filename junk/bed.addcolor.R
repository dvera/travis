bed.addcolor <-
function( bed, color, strand="." ){
	library(tools)
	ext<-file_ext(bed)
	curcol<-paste(col2rgb(color),collapse=",")
	curbed<-read.tsv(bed)
	outname<-paste0(removeext(bed),"_",color,".",ext)
	if(ncol(curbed)<3){ stop ("less than 3 columns in file")}
	if(ncol(curbed)==3){ curbed$V4 <- 1:nrow(curbed)}
	if(ncol(curbed)==4){ curbed$V5 <- 1}
	if(ncol(curbed)==5){ curbed$V6 <- strand}
	if(ncol(curbed)==6){ curbed$V7  <- curbed$V2}
	if(ncol(curbed)==7){ curbed$V8  <- curbed$V3}
	if(ncol(curbed)==8){ curbed$V9  <- curcol}
	write.tsv(curbed,file=outname)
}
