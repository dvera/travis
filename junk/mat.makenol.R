mat.makenol <-
function( bedfile, regionsize=3000, start=2, stop=3, suffix="", center=FALSE, row_names=TRUE, stepsize=1,path="/lustre/maize/home/dlv04c/data/b73v2/nol/a375/"){

	#REQUIRES NOL FILES TO BE SEPARATED BY CHROMOSOME, AND BE NAMED WITH SAME CHROMOSOME NAMES WITH .txt EXTENSION ( "chr1.txt" )
	
	pb <- txtProgressBar(min = 0, max = nrow(bedfile), style = 3)

	#REMOVE UNSHARED CHROMOSOMES FROM BED
	bedchroms<-unique(bedfile$V1)
	gffchroms<-unlist(lapply(strsplit(list.files(path=path),"\\."),"[",1))
	bed<-subset(bedfile,bedfile$V1 %in% gffchroms)
	
	#MINIMIZE REGION TABLE AND ADJUST CENTER FOR minus-stranded
	if(ncol(bed)>5){
		regions<-data.frame(V1=bed$V1,V2=bed[,start],V3=bed[,stop],V4=bed$V6,stringsAsFactors=FALSE)
		regions[which(regions[,4]=="-"),2]<-regions[which(regions[,4]=="-"),3]
	}
	else{
		if(center==TRUE){
		regions<-data.frame(V1=bed$V1,V2=round((bed$V2+bed$V3)/2),V3=bed$V3,V4=1)
		}
		else{
		regions<-data.frame(V1=bed$V1,V2=bed$V2,V3=bed$V3,V4=1)
		}
	}
	
	#CREATE MATRIX FOR STORING SCORES
	scorematrix<-matrix(ncol=regionsize,nrow=nrow(bed))
	if(row_names==TRUE){rownames(scorematrix)<-bed$V4}
	
	#FILL MATRIX WITH SCORES AND FLIP MINUS-STRANDED ROWS
	for(i in 1:nrow(regions)){
		startcoord<-round(regions[i,2] - regionsize/2)
		scorematrix[i,]<-system(command=paste("head -n ",startcoord+regionsize-25," ",path,regions[i,1],".txt | tail -n ",regionsize,sep=""),intern=TRUE)
		setTxtProgressBar(pb, i)
	}
	scorematrix[which(regions[,4]=="-"),1:regionsize]<-scorematrix[which(regions[,4]=="-"),regionsize:1]
	
	close(pb)
	
	#SAVE MATRIX FILE
	write.table(scorematrix,file=paste("nol-",suffix,".mat",collapse="",sep=""),sep="\t",quote=FALSE,row.names=TRUE,col.names=FALSE)
}
