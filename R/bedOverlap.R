bedOverlap <-
function( featurefiles,annotationfiles,suffix,genomefile,numshuffles=100,bpylim=5,b73=FALSE,genebed=NULL, scoremats=NULL, targetregions=NULL, cores="max", featnames=basename(removeext(featurefiles)), annonames=basename(removeext(annotationfiles)),usefeaturecenter=FALSE,useannocenter=FALSE,plotcolors=rainbow(length(featurefiles))){

	# setup multicore
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}

	# count number of input files
	numannos<-length(annotationfiles)
	numfeats<-length(featurefiles)

	# copy file names in new vectors
	feats<-featurefiles
	annos<-annotationfiles
	targets<-targetregions

	# check parameters
	if(length(featnames) != numfeats){stop("# of feature names must match # of features")}
	if(length(annonames) != numannos){stop("# of annotation names must match # of annotations")}
	if(is.null(suffix)){stop("must specify suffix")}

	# create a temporary directory for temp files
	tmpdir<-paste0("tmp_",paste(sample(letters,6),collapse=""),"/")
	dir.create(tmpdir)

	# count intervals
	cat("counting intervals\n")
	numtotalannos<-unlist(mclapply(annos,filelines,mc.cores=cores))
	numtotalfeats<-unlist(mclapply(feats,filelines,mc.cores=cores))

	# prune intervals to target regions
	cat("pruning annotations\n")
	if(is.null(targetregions) == FALSE){
		annos<-unlist(mclapply(annos,bedtoolsIntersect, b=targets, mc.cores=cores))
		feats<-unlist(mclapply(feats,bedtoolsIntersect, b=targets, mc.cores=cores))
	}

	# convert intervals to centers
	if(usefeaturecenter){
		cat("calculating feature centers\n")
		feats<-unlist(mclapply(feats,bedCenters,mc.cores=cores))
		sfeats<-unlist(mclapply(sfeats,bedCenters,mc.cores=cores))
	}
	if(useannocenter){
		cat("calculating annotation centers\n")
		annos<-unlist(mclapply(annos,bedCenters,mc.cores=cores))
	}

	# shuffle features
	cat("shuffling features\n")
	sfeatnames<-lapply(1:numfeats,function(f) unlist(lapply(1:numshuffles, function(s) paste0(tmpdir,featnames[f],"_shuffle_",s,".bed") ) ) )
	sfeats<-mclapply(1:numfeats,function(f) unlist(lapply(1:numshuffles,function(s) bedShuffle(featurefiles[f], genomefile=genomefile, outname=sfeatnames[[f]][s],include=if(is.null(targetregions)==FALSE){targetregions} else{NULL}))) , mc.cores=cores , mc.preschedule=F)

	# count intervals
	if(!is.null(targetregions)){
		cat("counting intervals\n")
		numeachanno<-unlist(mclapply(annos,filelines,mc.cores=cores))
		numeachfeat<-unlist(mclapply(feats,filelines,mc.cores=cores))
	} else{
		numeachanno<-numtotalannos
		numeachfeat<-numtotalfeats
	}






	# cat("finding nearby genes\n")
	# if(is.null(genebed)==FALSE){
	# 	feats<-unlist(mclapply(feats,bedClosest,bed2=genebed,mc.cores=cores))
	# 	sfeats<-unlist(mclapply(sfeats,bedClosest,bed2=genebed,mc.cores=cores))
	# }




	cat("finding overlaps\n")
	# find overlaps
	fxa<-as.data.frame(mclapply(1:numannos, function(a){
		unlist(lapply(1:numfeats,function(f){
			as.numeric(system(paste("bedtools intersect -u -a",feats[f],"-b",annos[a],"| wc -l"),intern=TRUE))
			}))
		},mc.cores=cores))

	# find shuffled overlaps
	sfxa<-as.data.frame(mclapply(1:numannos, function(a){
		unlist(lapply(1:numfeats,function(f){
			mean(as.numeric(unlist(lapply(1:numshuffles, function(z){
				system(paste("bedtools intersect -u -a",sfeats[[f]][z],"-b",annos[a],"| wc -l"),intern=TRUE)
			}))))
		}))
	},mc.cores=cores))

	# calculate summary statistics
	fxa.pct<-100*fxa/numeachfeat
	sfxa.pct<-100*sfxa/numeachfeat
	fxa.ratio<-log2(fxa.pct/sfxa.pct)
	fxa.ratio<-data.matrix(fxa.ratio)
	fxa.ratio[is.infinite(fxa.ratio)] <- NA


	# rename rows and columns
	row.names(fxa) = paste0("fxa_",featnames)
	row.names(sfxa) = paste0("sfxa_",featnames)
	row.names(fxa.pct) = paste0("fxa.pct_",featnames)
	row.names(sfxa.pct) = paste0("sfxa.pct_",featnames)
	rownames(fxa.ratio) = paste0("fxa.ratio_",featnames)
	colnames(fxa)=colnames(sfxa)=colnames(fxa.pct)=colnames(sfxa.pct)=colnames(fxa.ratio)=annonames

	# calculate label coordinates
	s<-rep(0,numfeats*numannos)
	s[(1:numfeats)*numannos+1]<-1
	s<-s[1:(numfeats*numannos)]

	s2<-rep(0,numfeats*numannos)
	s2[(1:numannos)*numfeats+1]<-1
	s2<-s2[1:(numfeats*numannos)]

	# start pdf connection
	pdf(file=paste("beddist_",suffix,".pdf",sep=""))

	# plot proportion of intervals overlapping with target regions
	par(mar=c(10,4,4,2))
	if(!is.null(targetregions)){
		barplot(c(numeachfeat,numeachanno)/c(numtotalfeats,numtotalannos),names=c(featnames,annonames),las=3,col=c(rep("black",numfeats),rep("grey50",numannos)), ylab="proportion of intervals overlapping with target regions")
	}

	# plot observed/expected enrichment
	par(mar=c(10,4,4,10),xpd=TRUE)
	barplot(unlist(as.data.frame(t(fxa.ratio))), col=rainbow(numannos),las=3,xaxt='n',space=s, ylab="feature enrichment log2(observed/shuffled)")
	axis(1,(1:numfeats)*ceiling(numannos+1)-ceiling(numannos/2),labels=featnames,las=3)
	legend(1+sum(s+1),max(as.vector(fxa.ratio)),legend=annonames,fill=rainbow(numannos))

	write.table(rbind(fxa,fxa.pct,sfxa.pct,fxa.ratio),file=paste0("beddist_",suffix,".tsv"),row.names=T,col.names=T,sep="\t",quote=F)
	# fxrat<-as.data.frame(apply(fxa.ratio,2,sort))
	# row.names(fxrat)=featnames
	# colnames(fxrat)=annonames
	#
	# barplot(unlist(as.data.frame((fxrat))), col=rainbow(numfeats),las=3,xaxt='n',space=s2)
	# axis(1,(1:numannos)*ceiling(numfeats+1)-ceiling(numfeats/2),labels=annonames,las=3)
	# legend(1+sum(s+1),max(as.vector(fxrat)),legend=featnames,fill=rainbow(numfeats))
	#
	# fxrat2<-t(as.data.frame(apply(as.data.frame(fxa.ratio),1,sort)))
	# row.names(fxrat2)=featnames
	# colnames(fxrat2)=annonames
	#
	# barplot(unlist(as.data.frame(t(fxrat2))), col=rainbow(numannos),las=3,xaxt='n',space=s, ylab="feature enrichment log2(observed/shuffled)")
	# axis(1,(1:numfeats)*ceiling(numannos+1)-ceiling(numannos/2),labels=featnames,las=3)
	# legend(1+sum(s+1),max(as.vector(fxrat2)),legend=annonames,fill=rainbow(numannos))



	barplot(unlist(as.data.frame((fxa.ratio))), col=rainbow(numfeats),las=3,xaxt='n',space=s2)
	axis(1,(1:numannos)*ceiling(numfeats+1)-ceiling(numfeats/2),labels=annonames,las=3)
	legend(1+sum(s+1),max(as.vector(fxa.ratio)),legend=featnames,fill=rainbow(numfeats))


	# barplot(unlist(as.data.frame(t(fxa.pct))), col=rainbow(numannos),las=3,xaxt='n')
	# axis(1,(1:numfeats)*numannos-numannos/2,labels=featnames)
	#
	# #barplot(unlist(as.data.frame(t(sfxa.pct))), names=rep(annonames,numfeats) , las=3)
	#
	# barplot(unlist(as.data.frame(t(fxa.ratio))), names=rep(annonames,numfeats) , las=3)
	#
	# a<-data.matrix(cbind(fxa.pct,sfxa.pct))
	# b<-data.matrix(fxa.pct/sfxa.pct)
	# b[is.infinite(b)]<-NA
	#
	# par(mar=c(10,4,4,2))
	# barplot(a,beside=T,las=2,ylab="% features overlapping with annotation",ylim=c(0,100),col=rainbow(numfeats))
	# legend("topright",legend=featnames,col=rainbow(numfeats),lwd=3)
	#
	# bp<-barplot(b,beside=T,las=2,ylab="enrichment over shuffled",col=rainbow(numannos))
	# bp<-barplot(b,beside=T,las=2,ylab="enrichment over shuffled",col=rainbow(numannos),ylim=c(0,bpylim))
	# legend("topright",legend=annonames,col=rainbow(numannos),lwd=3)
	# abline(h=1,col="grey30",lwd=1)
	# text(x=bp,y=b,labels=as.numeric(b),xpd=TRUE,srt=90,adj=0)
	# cat("check 1\n")


	# if(is.null(scoremats) == FALSE & is.null(genebed)==FALSE){
	# 	cat("pass1\n")
	# 	matlist<-mclapply(scoremats,read.mat,mc.cores=cores)
	# 	#if(do.call(identical,lapply(matlist,nrow))==FALSE){stop("matrices have different numbers of rows\n")}
	# 	if(b73==TRUE){
	# 		matgenes<-row.names(matlist[[1]])
	# 	} else{matgenes<-unlist(lapply(strsplit(row.names(matlist[[1]]),"_"),"[",1))}
	#
	# 	fxa.genes<-lapply(1:numfeats, function(f) {
	# 		lapply(1:numannos, function(a) {
	# 			unique(readLines(pipe(paste("cut -f 4",fxa[a,f]))))
	# 		})
	# 	})
	# 	sfxa.genes<-lapply(1:numfeats, function(f) {
	# 		lapply(1:numannos, function(a) {
	# 			unique(readLines(pipe(paste("cut -f 4",sfxa[a,f]))))
	# 		})
	# 	})
	# 	fnxa.genes<-lapply(1:numfeats, function(f) {
	# 		lapply(1:numannos, function(a) {
	# 			unique(readLines(pipe(paste("cut -f 4",fnxa[a,f]))))
	# 		})
	# 	})
	#
	# 	sfnxa.genes<-lapply(1:numfeats, function(f) {
	# 		lapply(1:numannos, function(a) {
	# 			unique(readLines(pipe(paste("cut -f 4",sfnxa[a,f]))))
	# 		})
	# 	})
	#
	#
	# 	all.genes<-unique(readLines(pipe(paste("cut -f 4",genebed))))
	# 	all.featgenes<-lapply(1:numfeats,function(x) unique(readLines(pipe(paste("cut -f 4",feats[x])))))
	# 	all.sfeatgenes<-lapply(1:numfeats,function(x) unique(readLines(pipe(paste("cut -f 4",sfeats[x])))))
	# 	cat("checking if scoremats has values\n")
	# 	for(s in 1:length(scoremats)){
	# 		cat("finding gene scores for scoremats",s,"\n")
	# 		fxa.genescores<-lapply(1:numfeats, function(f) {
	# 			lapply(1:numannos, function(a) {
	# 				if(fxa.counts[a,f]>0){
	# 					rowMeans(as.matrix(matlist[[s]][which(matgenes %in% fxa.genes[[f]][[a]]),]),na.rm=TRUE)
	# 				}
	# 				else{
	# 					NA
	# 				}
	# 			})
	# 		})
	# 		sfxa.genescores<-lapply(1:numfeats, function(f) {
	# 			lapply(1:numannos, function(a) {
	# 				if(sfxa.counts[a,f]>0){
	# 					rowMeans(as.matrix(matlist[[s]][which(matgenes %in% sfxa.genes[[f]][[a]]),]),na.rm=TRUE)
	# 				}
	# 				else{
	# 					NA
	# 				}
	# 			})
	# 		})
	# 		fnxa.genescores<-lapply(1:numfeats, function(f) {
	# 			lapply(1:numannos, function(a) {
	# 				if(fnxa.counts[a,f]>0){
	# 					rowMeans(as.matrix(matlist[[s]][which(matgenes %in% fnxa.genes[[f]][[a]]),]),na.rm=TRUE)
	# 				}
	# 				else{
	# 					NA
	# 				}
	# 			})
	# 		})
	#
	#
	# 		all.genescores<-rowMeans(matlist[[s]][which(matgenes %in% all.genes),],na.rm=TRUE)
	#
	# 		all.featgenescores<-lapply(1:numfeats,function(x) rowMeans(matlist[[s]][which(matgenes %in% all.featgenes[[x]]),],na.rm=TRUE))
	# 		all.sfeatgenescores<-lapply(1:numfeats,function(x) rowMeans(as.matrix(matlist[[s]][which(matgenes %in% all.sfeatgenes[[x]]),]),na.rm=TRUE))
	#
	#
	# 		#all<-list(1:numfeats, function(f) lapply(c(fxa.genescores,fnxa.genescores,sfxa.genescores), "[" ,  f ) )
	# 		#all<-c(list(all.genescores),all.featgenescores,unlist(fxa.genescores,recursive=F),unlist(fnxa.genescores,recursive=F))
	#
	# 		all<-unlist(lapply(1:numfeats,function(f) unlist(lapply(1:numannos,function(a) list(fxa.genescores[[f]][[a]],fnxa.genescores[[f]][[a]],sfxa.genescores[[f]][[a]])),recursive=F)),recursive=F)
	# 		#all<-unlist(lapply(1:numfeats,function(f) unlist(lapply(c(fxa.genescores,fnxa.genescores,sfxa.genescores), lapply , "[" ,  f ),recursive=F)),recursive=F)
	# 		all<-c(list(all.genescores),all.featgenescores,all.sfeatgenescores,all)
	# 		allnames<-basename(removeext(unlist(t(cbind(unlist(fxa),unlist(fnxa),unlist(sfxa))))))
	# 		ats=c(1,1+(1:(numfeats*2)),unlist(lapply(1:(numfeats*numannos),function(x) x*3+1:3+x))+2*numfeats-2)
	# 		par(mar=c(10,4,4,2))
	# 		b<-boxplot(all,plot=F)
	# 		boxplot(all,at=ats,col=c("white",rainbow(numfeats),rainbow(numfeats),rep(rainbow(numannos*numfeats),each=3)),names=c("all genes",featnames,featnames,allnames),las=3,cex.axis=0.7,outline=F,main=basename(removeext(scoremats[s])))
	# 		#axis,1,at=ats,labels=
	# 		points(ats,unlist(lapply(all,mean,na.rm=TRUE)))
	# 		ats2<-b$stats[5,]
	# 		text(ats,ats2,labels=paste(" ",lapply(all,length)),xpd=T,cex=0.7,srt=90,adj=0)
	#
	# 	}



	dev.off()
}
