scrapcode <-
function( ){
#if(mincovwin>0){
#				cat("finding poorly-covered features\n")
#				goodcovrows<-which(unlist(lapply(1:nrow(maskmat),function(x) length(which(maskmat[x,] > 0))))>mincovwin)
#				maskmat<-maskmat[goodcovrows,]
#				goodcovgenefile<-paste(basename(removeext(beds[i])),"_",basename(removeext(maskbed)),"_gc",mincovwin,"-",numwindows,".bed",sep="")
#				cat("saving well-covered regions in",goodcovgenefile,"\n")
#				write.tsv(bed[goodcovrows,],file=goodcovgenefile)
#
}
