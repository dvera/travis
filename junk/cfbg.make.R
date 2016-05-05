cfbg.make <-
function( bedfile ){
	bedname<-basename(removeext(bedfile))
	outname<-paste(bedname,".cfbg",sep="")
	cat(bedname,": calculating fragment centers\n")
	system(paste("awk '{a=int(($2+$3)/2+0.5); $4=$3-$2; $2=a; $3=a+1;print}' OFS='\t' ",bedfile," | sort -k1,1 -k2,2n > ",outname,sep=""))
	return(outname)
}
