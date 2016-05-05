anno.dist <-
function( ){
	#define intersect function
	
	
	
	
	
	xsect			<- function( feat,annos){
		write.tsv(feat[,1:3],file="tmpfeatures.tmp")
		annotations<-unique(annos$V4)
		numannos<-length(annotations)
		annocounts<-vector(length=numannos)
		for(i in 1:numannos){
			curanno<-subset(annos,annos$V4==annotations[i])
			write.tsv(curanno[,1:5],file="tmpanno.tmp")
			sc("bedtools intersect -u -a tmpfeatures.tmp -b tmpanno.tmp > tmpintersect.tmp")
			cat("pass3\n")
			if(file.info("tmpintersect.tmp")$size!=0){
				annocounts[i]<-nrow(read.tsv("tmpintersect.tmp"))
			}
			else{
				annocounts[i]<-0
			}
			cat("pass4\n")
		}
		annocounts
	}
	

	
	
	#prepare annotation file
	#allow no results to exist to complete script
	if(annotation.type=="gff"){
		cat("\tConverting gff to pga\n")
		annotation<-make.pga(annotation,intern=TRUE)
	}
	if(annotation.type=="pga"){
	}
	#setup results table
	cat("\n\tPreparing data\n")
	annotations<-unique(annotation$V4)
	numannos<-length(annotations)
	results<-matrix(nrow=length(featurelist),ncol=numannos)
	expresults<-matrix(nrow=length(featurelist),ncol=numannos)
	colnames(results)<-annotations
	row.names(results)<-removeext(featurelist)
	resultsd<-results
	annoave<-results
	
	
	write.tsv(annotation,file="tmpannotation.tmp")
	expmat<-as.data.frame(read.mat(expressionmat)[,1],stringsAsFactors=FALSE)
	expresults<-matrix(ncol=numannos,nrow=nrow(expmat))

	
	
	
	
	
	#intersect each feature with annotations
	cat("\n\tIntersecting features with annotation")
	write.tsv(annotation,file="tmpannos.tmp")
	filler<-as.data.frame(matrix(nrow=nrow(annotation),ncol=length(featurelist)))
	annotation2<-cbind(annotation[,4],annotation[,5],filler)
	colnames(annotation2)[1:2]<-c("type","gene")
	for(f in 1:length(featurelist)){
		
		cat("\n\t\t",featurelist[f])
		features<-read.tsv(featurelist[f])
		if(feature.center==TRUE){
			features$V2<-floor((features$V3+features$V2)/2)
			features$V3<-features$V2+1
		}
		cat("pass\n")
		#get average feature scores for each annotation
		results[f,]<-xsect(features,annotation)
		cat("pass2\n")
		#get average expression for each annotation overlapping with each feature
		sc(paste("bedtools intersect -u -a tmpannotation.tmp -b ",featurelist[f], " > tmpintersect2.tmp",sep=""))
		annotation3<-read.tsv("tmpintersect2.tmp")
		annotation3<-annotation3[,4:5]
		annotation3<-unique(annotation3)
		annotation4<-merge(annotation3,expmat,by.x="V5",by.y="row.names")
		for(i in 1:numannos){
			curanno<-subset(annotation4,annotation4[,2]==annotations[i])
			annoave[f,i]<-mean(curanno[,3],na.rm=TRUE)
		}
	}
	
	#determine gene expression based on gene annotation overlap with features
	sc(paste("bedtools intersect -c -a tmpannos.tmp -b ",featurelist[f], " > tmpintersect3.tmp",sep=""))

	
	#print(head(annotation2,n=50))
	for(i in 1:numannos){
		annoscores<-subset(annotation2,annotation2[,2]==annotations[i])
		write.tsv(annoscores,file="tmpannos2.tmp")
		scoregenes<-unique(annoscores[,3])
		numgenes<-length(scoregenes)
	}
	
	
	cat("\n\t\t Denominator")
	if(length(denominator)>1){
		#denom<-as.vector(read.table(paste(denominator),header=TRUE,row.names=1,stringsAsFactors=FALSE))
		denom<-denominator
		for(i in 1:length(featurelist)){
			resultsd[i,]<-(results[i,]/denom)*1000
		}
	}
	
	
	if(length(denominator)==1){
		features<-read.tsv(denominator)
		if(feature.center==TRUE){
			features$V2<-floor((features$V3+features$V2)/2)
			features$V3<-features$V2+1
		}
		write.tsv(features[,1:3],file="tmpfeatures.tmp")
		denom<-xsect(features,annotation)
		write.table(denom,file="denominator.tsv",sep="\t",quote=FALSE)
		for(i in 1:length(featurelist)){
			resultsd[i,]<-(results[i,]/denom)*1000
		}
	}
	
	
	
	cat("\n")
	write.table(results,file=paste(prefix,"-annodist.txt",sep=""),quote=FALSE,sep="\t")
	
	if(X==FALSE){pdf(file=paste(prefix,"-annodist.pdf",sep=""))}
	par(mfrow=c(3,3))
	barplot(results, beside=TRUE,las=3,col=plotcolors,axis.lty=0,xlab="Region",ylab="# Features")
	legend("topleft",legend=legendnames, fill=plotcolors)

	barplot(resultsd, beside=TRUE,las=3,col=plotcolors,axis.lty=0,xlab="Region",ylab="Features / kb")
	legend("topleft",legend=legendnames, fill=plotcolors)

	barplot(annoave, beside=TRUE,las=3,col=plotcolors,axis.lty=0,xlab="Region",ylab="Expression Level")
	legend("topleft",legend=legendnames, fill=plotcolors)

	
	results<-results[,colnames(results) %in% pieFeatures]
	resultsd<-resultsd[,colnames(resultsd) %in% pieFeatures]
	for(y in 1:nrow(results)){
		pie(results[y,],labels=paste(colnames(results),",",results[y,]),clockwise=TRUE,main=paste("Total",row.names(results)[y]))
	}
	plot(0,0)
	for(y in 1:nrow(resultsd)){
		pie(resultsd[y,],labels=paste(colnames(resultsd),",",resultsd[y,]),clockwise=TRUE,main=paste("Features/kb",row.names(results)[y]))
	}
	if(X==FALSE){dev.off()}
}
