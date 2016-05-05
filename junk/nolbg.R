nolbg <-
function( bg ){
	#assumes no data within 25 bp of end of chromosome
	#assumes nol file names match column1 of bedgraph and have .txt extension
	chroms<-unique(bg$V1)
	bg$V1<-bg$V1+25
	for(i in 1:length(chroms)){
		write.tsv(bg[bg$V1==chroms[i],][,1],file="linenumbers.tmp")
		sc(paste("awk 'NR==FNR { a[$1];next } (FNR in a)' linenumbers.tmp ",chroms[i],".txt > ",chroms[i],".tmp",sep=""))
		
	}
}
