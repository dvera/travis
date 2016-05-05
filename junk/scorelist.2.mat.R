scorelist.2.mat <-
function( scorefiles, genebed, namecol=1, delimiter="", namenumber=1, scorecol=4, nametype="symbol",matcol=200){
	
	bed<-read.tsv(genebed)
	if(nametype=="symbol"){bedids<-bed$V13}
	if(nametype=="ucsc"){bedids<-bed$V4}
	bedids<-data.frame("V1"=bedids,stringsAsFactors=FALSE)
	matrownames<-paste(bed$V4,bed$V1,bed$V13,sep="_")
	numscores<-length(scorefiles)

	scores<-lapply(scorefiles,read.tsv)
	scores<-lapply(1:numscores, function(s){
		if(delimiter!=""){
			scores[[s]][,namecol]<-unlist(lapply(strsplit(scores[[s]][,namecol],delimiter, fixed=TRUE),"[",namenumber))
		}
		scores[[s]]<-scores[[s]][-which(duplicated(scores[[s]][,namecol])),]
		scoretable<-data.frame("V1"=scores[[s]][,namecol],"V2"=scores[[s]][,scorecol],stringsAsFactors=FALSE)
		scoretable
	})
	
	scoremats<-lapply(1:numscores, function(s){
		scoremat<-merge(bedids,scores[[s]],by="V1",all.x=TRUE,all.y=FALSE,sort=FALSE)
		scoremat<-scoremat[match(scoremat$V1,bedids$V1),]
		mat<-matrix(NA_real_,nrow=nrow(bed),ncol=matcol)
		mat[,]<-scoremat$V2
		row.names(mat)<-matrownames
		mat
	})
	
	unlist(lapply(1:numscores, function(s){
		write.mat(scoremats,file=paste(basename(removeext(scorefiles[s])), "_", basename(removeext(genebed)),".mat10", sep=""))
	}))
}
