bedRescore <- #rename to bedRescore
function ( scores , bedfile , countCols , nameCol=1 , samplenames=NULL ){

	options(scipen=9999)
	bed=tsvRead(bedfile)
	cnt=tsvRead(scores,col_names=T,comment="#")
	bedname <- basename(removeext(bedfile))
	bedcols<-ncol(bed)
	bedrows<-nrow(bed)
	#matgenes<-unlist(lapply(strsplit(row.names(mat),";") , "[" , 3 ))

	# if(b73){
	# 	matgenes<-remove.suffix(matgenes,"_T")
	# 	matgenes<-gsub("_FGT","_FG",matgenes)
	# }
	#

	namematch<-match(bed[,4],cnt[,nameCol])
	if(length(na.omit(namematch))==0){stop("no name matches between bed and scores")}
	cnt <- cnt[namematch,]

	numsamples<-length(countCols)
	samplecols<-countCols

	if(is.null(samplenames)){ samplenames<-basename(colnames(cnt)[countCols]) }

	outnames <- paste0(bedname,"_",samplenames,".bed")

	dump <- lapply(seq_len(numsamples),function(x){
		curbed <- bed
		curbed[,5] <- as.numeric(cnt[,samplecols[x]])
		#rpkms<-rpkm.default(readcounts,cnt[,6])

		#mat <- matrix(rpkms,nrow=nrow(bed),ncol=1)
		#row.names(mat) <- matrownames
		tsvWrite(curbed[order(curbed[,1],curbed[,2],curbed[,3]),], outnames[x])
	})

	return(outnames)



}
