bedWords <-
function( bedfiles, genomefa , pdfname , numbases=50, sizerange=c(1,Inf), numfrags=NULL, reference="center", symmetric=TRUE , cores="max", lspan=0, slop=0 , strand=FALSE, plotcolors=NA , legendnames=NA, acgt.colors=c("red","purple","blue","green") ){

	if(grepl(".pdf",pdfname)==FALSE){stop("pdf name must be defined with .pdf extension")}
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}

	if(reference %in% c("center","end") == FALSE){stop("not a valid reference point, use 'center' or 'end'")}
	if(is.na(legendnames)){ legendnames <- basename(removeext(bedfiles)) }
	tmpdir <- paste0(basename(tempdir()),"/")
	dir.create (tmpdir)

	numbeds<-length(bedfiles)
	beds<-bedfiles
	bednames<-legendnames


	if(is.na(plotcolors)==FALSE & length(plotcolors) != numbeds){ stop("lengths of plotcolors and number of beds should be equal")}
	if(cores > numbeds){ cores <- numbeds }

	cat("determining dinucleotide words\n")
	nwords<-expand.grid(c("A","C","G","T"),c("A","C","G","T"))
	nnums<-expand.grid(1:4,1:4)
	dwords<-paste0(nwords[,1],nwords[,2])
	numwords<-length(dwords)


	beds<-unlist(mclapply(1:numbeds, function(i){


		if(is.infinite(sizerange[2]) == FALSE & sizerange[1]>1){
			beds[i]<-bedParseLengths(beds[i],brks=sizerange)
		}

		if(is.null(numfrags) == FALSE){
			cat(bednames[i],": grabbing",numfrags,"random fragments\n")
			newname<-paste(tmpdir,bednames[i],"_",numfrags,"random.bed",sep="")
			system(paste("shuf -n",numfrags,beds[i],">",newname))
			beds[i]<-newname
		}

		if(reference=="center"){
			cat(bednames[i],": calculating fragment centers\n")
			newname<-paste(tmpdir,bednames[i],"_centersforwords.bed",sep="")
			system(paste("awk 'BEGIN{b=(",numbases,"/2)+slop} {a=int(($2+$3)/2)-b; $2=a; $3=a+",numbases,"+slop;print}' OFS='\t' ",beds[i],">",newname))
			beds[i]<-newname
		} else{
			cat(bednames[i],": calculating fragment ends\n")
			newname<-paste(tmpdir,bednames[i],"_endsforwords.bed",sep="")
			system(paste("awk '{$2=$2-",slop,"; $3=$2+",numbases,"+slop;print}' OFS='\t' ",beds[i],">",newname))
			beds[i]<-newname
		}

		return(beds[i])
	},mc.cores=detectCores()))


	print(beds)

	basemats<-mclapply(1:numbeds, function(i){
		cat(bednames[i],": extracting sequences from bed\n")
		seqbedname<-paste0(tmpdir,bednames[i],".seqbed")
		basemat<-toupper(readLines(pipe(paste("bedtools getfasta -tab",if(strand){"-s"},"-fi",genomefa,"-bed",beds[i],"-fo stdout | cut -f 2" ) ) ) )

		print(head(basemat))
		print(length(basemat))

		numseqs<-length(basemat)
		basemat<-t(simplify2array(lapply(1:numseqs,function(x) substring(basemat[x],1:numbases,1:numbases))))
		basemat[basemat=="N"]<-NA
		as.data.frame(basemat)
	},mc.cores=cores)

	#basemats<-lapply(basemats,as.data.frame)
	# basemats2<-as.matrix(basemats[[1]])
	# basemats2[basemats2=="A"]<--1
	# basemats2[basemats2=="T"]<--1
	# basemats2[basemats2=="C"]<-1
	# basemats2[basemats2=="G"]<-1
	# mode(basemats2)<- "numeric"
	# image(t(basemats2[order(kmeans((basemats2),2)$cluster),]))


	nucfreqs<-mclapply(1:numbeds, function(i){
		freqs<-as.data.frame(lapply(1:ncol(basemats[[i]]),function(x) as.vector(table(basemats[[i]][,x]))))
		freqsums<-matrix(rep(apply(freqs,2,sum),4),nrow=4)
		freqs<-freqs/freqsums
		row.names(freqs)<-c("A","C","G","T")
		colnames(freqs)<-1:ncol(freqs)
		freqs
	},mc.cores=cores)




	expectedwords<-lapply(1:numbeds, function(b){

		ew<-as.data.frame(mclapply(1:(numbases-1),function(n){
			unlist(lapply(1:16,function(x){
				nucfreqs[[b]][nnums[x,1],n]*nucfreqs[[b]][nnums[x,2],n+1]
			}))
		},mc.cores=cores))
		colnames(ew)<-1:(numbases-1)

		ew<-rbind(ew,colSums(ew[c(1,4,13,16),]))
		ew<-rbind(ew,colSums(ew[c(6,7,10,11),]))

		row.names(ew)<-c(dwords,"TT/AA/TA/AT","CC/GG/CG/GC")

		colnames(ew)<-1:(numbases-1)

		ew
	})




	fragnums<-unlist(lapply(basemats,nrow))

	freqmats<-list()
	for(j in 1:numbeds){

		cat(bednames[j],": finding dinucleotide frequencies\n")
		freqmat<-simplify2array(mclapply(1:(numbases-1), function(x) {
			curwords<-paste(basemats[[j]][,x], basemats[[j]][,x+1], sep="")
			unlist(lapply(1:numwords, function(y) length(which(curwords==dwords[y]))))
		},mc.cores=detectCores()))
		numseqs<-nrow(basemats[[j]])
		#freqmat<-(freqmat/nrow(basemats[[j]]))*100
		freqmat<-rbind(freqmat,colSums(freqmat[c(1,4,13,16),]))
		freqmat<-rbind(freqmat,colSums(freqmat[c(6,7,10,11),]))
		row.names(freqmat)<-c(dwords,"TT/AA/TA/AT","CC/GG/CG/GC")
		freqmats[[j]]<-freqmat
		#cat("saving dinucleotide frequency matrix\n")
		#write.mat(freqmat,file=paste(bednames[j],".dnfmat",sep=""))

	}

	freqmats<-lapply(1:numbeds,function(x) (freqmats[[x]]/fragnums[x]))
	adjfreqmats<-lapply(1:numbeds,function(x) log2(freqmats[[x]]/expectedwords[[x]]))

	if(symmetric){
		freqmats<-lapply(1:numbeds,function(x) (freqmats[[x]][,1:(numbases-1)] + freqmats[[x]][,(numbases-1):1])/2)
		adjfreqmats<-lapply(1:numbeds,function(x) (adjfreqmats[[x]][,1:(numbases-1)] + adjfreqmats[[x]][,(numbases-1):1])/2)
		nucfreqs<-lapply(1:numbeds,function(x) (nucfreqs[[x]][,1:(numbases)] + nucfreqs[[x]][,(numbases):1])/2)
	}

	nucfreqs<-lapply(nucfreqs,t)
	nucfreqs<-lapply(nucfreqs,as.data.frame)
	freqmats<-lapply(freqmats,t)
	freqmats<-lapply(freqmats,as.data.frame)
	adjfreqmats<-lapply(adjfreqmats,t)
	adjfreqmats<-lapply(adjfreqmats,as.data.frame)

	pdf(file=pdfname,width=10,height=10)

	if(is.na(plotcolors) ) { rb<-rainbow(numbeds) }


	cat("plotting absolute nucleotide frequencies\n")
	par(mfrow=c(5,1),mar=c(0,3,0,1),oma=c(4,4,4,4),xpd=TRUE)
	plot(0,type="n", axes=FALSE , ann=FALSE)
	legend("topleft",legend=paste(bednames,"(",fragnums," fragments)"),col=rb,lwd=2)
	for(j in 1:4){
		all<-unlist(lapply(nucfreqs,"[",j))
		plot(0,type="n",xlim=c(0,numbases)-slop,ylim=c(min(all,na.rm=TRUE),max(all*1.2,na.rm=TRUE)),xaxt='n',ylab="proportion of all dinucleotides",xlab="nucleotide")
		abline(v=seq(from=0,to=numbases-slop,by=10),lwd=1,col="grey80")

		# make vertical lines at references
		if(reference=="center"){
			abline(v=c(0,numbases/2),lwd=1,col="black")
		} else{
			abline(v=0,lwd=1,col="black")
		}


		for(k in 1:numbeds){
			x<-(1:(numbases))-slop
			y<-unlist(nucfreqs[[k]][,j])
			if(lspan > 0){
				y<-loess(y~x,span=lspan)$fitted
				if(symmetric){
					y<-(y[1:length(y)]+y[length(y):1])/2
				}
			}
			lines(x,y,col=rb[k])
		}


		if(j==4){axis(side=1,at=seq(from=-slop,to=numbases-slop,by=10))}
		#legend("topright",legend=colnames(nucfreqs[[1]][j]))
		mtext(colnames(nucfreqs[[1]][j]),side=4,line=1,las=1,cex=2)

	}
	mtext("absolute nucleotide frequencies (proportion)",outer=T,side=2)

	par(mfrow=c(4,1),mar=c(0,3,1,4),oma=c(4,4,4,4),xpd=TRUE)

	for(j in 1:numbeds){
		plot(0,type="n",xlim=c(-slop,numbases-slop),ylim=c(min(unlist(nucfreqs[[j]])),max(unlist((nucfreqs[[j]])))),main=bednames[j],xaxt='n',ylab="log2(observed/expected)",xlab="nucleotide")
		abline(v=seq(from=-slop,to=numbases-slop,by=10),lwd=1,col="grey80")

		# make vertical lines at references
		if(reference=="center"){
			abline(v=c(0,numbases/2),lwd=1,col="black")
		} else{
			abline(v=0,lwd=1,col="black")
		}
		x<-((1:numbases))-slop
		lines(x,nucfreqs[[j]][,1],type="l",col=acgt.colors[1])
		lines(x,nucfreqs[[j]][,2],type="l",col=acgt.colors[2])
		lines(x,nucfreqs[[j]][,3],type="l",col=acgt.colors[3])
		lines(x,nucfreqs[[j]][,4],type="l",col=acgt.colors[4])
		legend("topright",inset=c(-0.05,0),legend=colnames(nucfreqs[[1]]),col=acgt.colors,lwd=3)
		if(j==numbeds){axis(side=1,at=seq(from=-slop,to=numbases-slop,by=10))}

		#mtext(bednames[j],side=4)
	}
	mtext("absolute nucleotide frequencies (proportions)",outer=T,side=2)

	cat("plotting absolute dinucleotide frequencies\n")
	par(mfrow=c(4,4),mar=c(0,1,0,1),oma=c(4,4,4,4),xpd=T)
	for(j in 1:18){
		all<-unlist(lapply(freqmats,"[",j))
		plot(0,type="n",xlim=c(-slop,numbases-slop),ylim=c(min(all,na.rm=TRUE),max(all,na.rm=TRUE)),xaxt='n',ylab="")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")

		if(reference=="center"){
			abline(v=c(0,slop+numbases/2),lwd=1,col="black")
		} else{
			abline(v=0,lwd=1,col="black")
		}

		for(k in 1:numbeds){
			x<-(1:(numbases-1))-slop
			y<-freqmats[[k]][,j]
			if(lspan > 0){
				y<-loess(y~x,span=lspan)$fitted
				if(symmetric){
					y<-(y[1:length(y)]+y[length(y):1])/2
				}
			}
			#lines(x,y,type="s",col=rgb(0,0,0,75,maxColorValue=255)
			lines(x,y,col=rb[k])
		}
		#legend("topleft",legend=paste(bednames,"(",fragnums," fragments)"),col=rb,lwd=2)
		legend("top",legend=row.names(freqmat)[j],bty='n')
		if(j %in% 13:18){ axis(side=1,at=seq(from=-slop,to=numbases-slop,by=10)) }
		if(j == 16){
			mtext("absolute dinucleotide frequencies (proportions)",outer=T)
			par(mfrow=c(2,1),mar=c(2,4,2,4),oma=c(4,4,4,4),xpd=F)
		}
	}
	mtext("absolute dinucleotide frequencies (proportions)",outer=T)


	cat("plotting adjusted dinucleotide frequencies\n")
	par(mfrow=c(4,4),mar=c(0,1,0,1),oma=c(4,4,4,4))
	for(j in 1:18){
		all<-unlist(lapply(adjfreqmats,"[",j))
		all[is.infinite(all)] <- NA
		ylims<-c(min(all,na.rm=TRUE),max(all,na.rm=TRUE))
		print(ylims)
		plot(0,type="n",xlim=c(-slop,numbases-slop),ylim=ylims,xaxt='n',ylab="")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")

		if(reference=="center"){
			abline(v=c(slop,slop+numbases/2),lwd=1,col="black")
		} else{
			abline(v=slop,lwd=1,col="black")
		}

		for(k in 1:numbeds){
			x<-(1:(numbases-1))-slop
			y<-adjfreqmats[[k]][,j]
			y[is.infinite(y)] <- NA
			if(lspan > 0){
				y<-loess(y~x,span=lspan)$fitted
				if(symmetric){
					y<-(y[1:length(y)]+y[length(y):1])/2
				}
			}
			#lines(x,y,type="s",col=rgb(0,0,0,75,maxColorValue=255)
			lines(x,y,col=rb[k])
		}
		#legend("topright",legend=paste(bednames,"(",fragnums," fragments)"),col=rb,lwd=2)
		legend("top",legend=row.names(freqmat)[j],bty='n')
		if(j %in% 13:18){ axis(side=1,at=seq(from=-slop,to=numbases-slop,by=10)) }
		if(j == 16){
			mtext("observed/expected nucleotide frequencies (proportions)",outer=T)
			par(mfrow=c(2,1),mar=c(2,4,2,4),oma=c(4,4,4,4))
		}
	}
	mtext("observed/expected nucleotide frequencies (proportions)",outer=T)
	dev.off()
}
