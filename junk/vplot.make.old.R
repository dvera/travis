vplot.make.old <-
function( fragments,bed,regionsize=1000, suffix="",ylims=c(0,200),hm.max="default",X=TRUE,savepng=TRUE){

	#get object names for file names
	f<-deparse(substitute(fragments))
	b<-deparse(substitute(bed))

	#define distances from annotations to search for features
	leftboundary<--0.5*regionsize
	rightboundary<-0.5*regionsize

	#load necessary libraries
	library(gplots)
	library(compiler)
	
	#define feature-mapping function
	f1 <- function(frags, centers) {
		score_matrix <- matrix(NA_real_, 1, 2)
		d1 <- frags[,1]
		d2 <- frags[,2]
		for (i in seq_along(centers)) {
			score <- d1 - centers[i]
			idx <- score > leftboundary & score <= rightboundary #rows which have relative coords of -500 and 500
			if(length(idx)>0){
				score_matrix<-rbind(score_matrix,cbind(score[idx],d2[idx]))
			}
			setTxtProgressBar(pb, y)
			y=y+1
		}
		
		cat("\n")
		score_matrix
	}
	
	#compile f1
	f1c <- cmpfun(f1)

	
	cat("removing unshared chromosomes\n")
	fragments<-subset(fragments,fragments$V1 %in% bed$V1)
	bed<-subset(bed,bed$V1 %in% fragments$V1)

	
	cat("preprocessing data\n")
	#calculate fragment sizes
	fragments$V4=fragments$V3-fragments$V2
	#calculate center of fragments
	fragments$V2=floor( (fragments$V2+fragments$V3)/2 )
	#calculate center of annotations
	bed$V2=floor( (bed$V2+bed$V3)/2 )

	y=1
	chromosomes<-unique(bed[,1])
	
	cat("searching ",nrow(fragments)," fragments for those within ",regionsize/2," bp of ",nrow(bed)," annotations\n",sep="") 
	pb <- txtProgressBar(min = 0, max = nrow(bed), style = 3)
	for(c in 1:length(chromosomes)){
		#find fragments in chromosome
		chromfrags<-subset(fragments,fragments[,1]==chromosomes[c])
		#trim to relevant data
		chromfrags<-data.matrix(cbind(chromfrags[,2],chromfrags[,4]))
		#find annotations in chromosome
		chromregions<-subset(bed,bed[,1]==chromosomes[c])
		centers<-chromregions[,2]
		#find fragments nearby annotations
		chrommatrix<-f1c(chromfrags,centers)
		#create scorematrix if first chromosome processed
		if(c==1){scorematrix<-chrommatrix[2:nrow(chrommatrix),]}
		#append chromscores to scorematrix if not first chromosome processed
		else{scorematrix<-rbind(scorematrix,chrommatrix[2:nrow(chrommatrix),])}
	}
	close(pb)
	
	#define matrix for heatmap
	m<-matrix(0,ncol=regionsize,nrow=ylims[2])
	
	#remove too-distant fragments and trim regions extending beyond matrix
	scorematrix[,1]=scorematrix[,1]+regionsize/2
	scorematrix<-scorematrix[which(scorematrix[,1]<regionsize),]
	scorematrix<-scorematrix[which(scorematrix[,1]>0),]
	scorematrix<-scorematrix[which(scorematrix[,2]<=ylims[2]),]
	cat(nrow(scorematrix),"/",nrow(fragments)," within size and range of annotations\n")
	
	#fill matrix with data
	for(d in 1:nrow(scorematrix)){
		mrow=scorematrix[d,2]
		mcol=scorematrix[d,1]
		m[mrow,mcol]=m[mrow,mcol]+1
	}
	
	#flip matrix vertically
	m<-m[nrow(m):1,]
	
	#define maximum in heatmap
	if(hm.max=="default"){hm.max=as.numeric(round(quantile(m[m!=0],prob=0.99)))}
	
	#define breaks in heatmap colors
	brks<-seq(0,hm.max,by=hm.max/100)
	
	#define heatmap colors
	greycolors<-grey((brks   [1:(   length(brks) -1)]  )/hm.max)
	
	#define file name based on objects and parameters
	filename<-paste("vplot-",f,"-",b,"-r",regionsize,"-f",ylims[2],"-d",hm.max,suffix,collapse="",sep="")
	
	#save matrix for later replotting
	write.mat(m,file=paste(filename,".mat",sep=""))

	#open png
	if(savepng==TRUE){
		png(file=paste(filename,".png",sep=""),width=1000,height=1000)
	}
	
	#draw heatmap
	heatmap.2(m,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=seq(0,hm.max,by=hm.max/100),labRow=NA,labCol=NA,col=greycolors)
	
	#close png
	if(savepng==TRUE){
		dev.off()
	}

	if(X==TRUE){
		#draw heatmap
		heatmap.2(m,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=seq(0,hm.max,by=hm.max/100),labRow=NA,labCol=NA,col=greycolors)
	}
}
