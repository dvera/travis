vplot.make.old2 <-
function( fragments,bed,regionsize=1000, suffix="",ylims=c(0,200),hm.max="default"){
	leftboundary<--0.5*regionsize-2-ylims[2]/2
	rightboundary<-0.5*regionsize+2+ylims[2]/2
	library(gplots)
	library(compiler)
	f1 <- function(frags, centers) {
		score_matrix <- matrix(NA_real_, 1, 2)
		d1 <- frags[,1]
		d2 <- frags[,2]
		
		for (i in seq_along(centers)) {
			score <- d1 - centers[i]
			idx <- score > leftboundary & score <= rightboundary
			if(length(idx)>0){
				score_matrix<-rbind(score_matrix,cbind(score[idx],d2[idx]))
			}
			setTxtProgressBar(pb, i)
		}
		cat("\n")
		score_matrix
	}
	f1c <- cmpfun(f1)

	cat("removing unshared chromosomes\n")
	fragments<-subset(fragments,fragments$V1 %in% bed$V1)
	bed<-subset(bed,bed$V1 %in% fragments$V1)

	#create bed matrix with relevant data
	fragments$V4=fragments$V3-fragments$V2
	fragments$V2=floor( (fragments$V2+fragments$V3)/2 )
	bed$V2=floor( (bed$V2+bed$V3)/2 )

	cat("#regions=",nrow(bed),"\n")
	chromosomes<-unique(bed[,1])
	for(c in 1:length(chromosomes)){
		cat("chromosome",chromosomes[c],"\n")
		chromfrags<-subset(fragments,fragments[,1]==chromosomes[c])
		chromfrags<-data.matrix(cbind(chromfrags[,2],chromfrags[,4]))
		chromregions<-subset(bed,bed[,1]==chromosomes[c])
		centers<-chromregions[,2]
		pb <- txtProgressBar(min = 0, max = length(centers), style = 3)
		chrommatrix<-f1c(chromfrags,centers)
		if(c==1){scorematrix<-chrommatrix[2:nrow(chrommatrix),];cat("created scorematrix\n")}
		else{scorematrix<-rbind(scorematrix,chrommatrix[2:nrow(chrommatrix),])}
		cat("appended scorematrix\n")
	}
	
	
	
	
	#shift coordinates so that they are all positive for matrix
	scorematrix[,1]=scorematrix[,1]+regionsize/2
	scorematrix<-scorematrix[which(scorematrix[,2]<=ylims[2]),]
	
	#define left and right boundary distance for each fragment
	flanks<-round(scorematrix[,2]/2)
	
	#add fragment left and right boundaries
	scorematrix<-cbind(scorematrix,scorematrix[,1]-flanks,scorematrix[,1]+flanks)
	
	#remove fragments beyond heatmap boundaries
	scorematrix<-scorematrix[which(scorematrix[,3]<regionsize),]
	scorematrix<-scorematrix[which(scorematrix[,4]>0),]
	
	#truncate fragments extending past heatmap boundaries
	scorematrix[,3][scorematrix[,3]<0]=0
	scorematrix[,4][scorematrix[,4]>regionsize]=regionsize
	
	
	#print(nrow(scorematrix))
	#plot(density(scorematrix[,2],na.rm=TRUE))
	#X11()
		#plot(scorematrix[,1],scorematrix[,2],pch=20,col=rgb(0,0,0,5,maxColorValue=255))
		#X11()

	
	write.tsv(scorematrix,file="scoremat.tsv")
	for(d in 1:nrow(scorematrix)){
		mrow=scorematrix[d,2]
		mcol=scorematrix[d,1]
		m[mrow,(scorematrix[d,3]:scorematrix[d,4])]=m[mrow,(scorematrix[d,3]:scorematrix[d,4])]+1
	}
	
	#flip matrix vertically
	m<-m[nrow(m):1,]	

	
	if(hm.max=="default"){hm.max=as.numeric(round(quantile(m[m!=0],prob=0.99)))}
	brks<-seq(0,hm.max,by=hm.max/100)
	greycolors<-grey((brks   [1:(   length(brks) -1)]  )/hm.max)

	
	
	filename<-paste("vplot-",paste(deparse(substitute(fragments))),"-",paste(deparse(substitute(bed))),"-r",regionsize,"-f",ylims[2],"-d",hm.max,suffix,".mat",collapse="",sep="")
	write.mat(m,file=paste(filename,".mat",sep=""))
	if(png==TRUE){
		png(file=paste(filename,".png",sep=""))
	}
	heatmap.2(m,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=seq(0,hm.max,by=hm.max/100),labRow=NA,labCol=NA,col=greycolors)
	if(png==TRUE){
		dev.off()
	}
}
