bed.addstrand <-
function( bed, strand="+" ){
	ext<-file_ext(bed)
	curbed<-read.tsv(bed)
	if(strand=="+"){ strandname=="plus" }
	if(strand=="-"){ strandname=="minus"}
	outname<-paste0(removeext(bed),"_",strandname,".",ext)
	if(ncol(curbed)<3){ stop ("less than 3 columns in file")}
	if(ncol(curbed)==3){ curbed$V4 <- 1:nrow(curbed)}
	if(ncol(curbed)==4){ curbed$V5 <- 1}
	if(ncol(curbed)==5){ curbed$V6 <- strand}
	write.tsv(curbed,file=outname)
}
