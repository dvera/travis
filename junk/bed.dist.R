bed.dist <-
function( featurefiles,annotationfiles,suffix,numshuffles=100,bpylim=5,b73=FALSE,genebed=NULL, scoremats=NULL, targetregions=NULL, cores="max", featnames=basename(removeext(featurefiles)), annonames=basename(removeext(annotationfiles)),usefeaturecenter=FALSE,useannocenter=FALSE,plotcolors=rainbow(length(featurefiles))){
	
	#gather info
	numannos<-length(annotationfiles)
	numfeats<-length(featurefiles)
	feats<-featurefiles
	annos<-annotationfiles
	if(length(featnames) != numfeats){stop("# of feature names must match # of features")}
	if(length(annonames) != numannos){stop("# of annotation names must match # of annotations")}
	
	tmpdir<-paste0("tmp_",paste(sample(letters,6),collapse=""))
	dir.create(tmpdir)
	
	if(is.null(suffix)){stop("must specify suffix")}
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	if(is.null(targetregions) == FALSE){
		cat("merging target regions\n")
		targets<-bedtoolsMerge(targetregions)
		cat("counting total annotations\n")
		numtotalannos<-unlist(mclapply(annos,filelines,mc.cores=cores))
		numtotalfeats<-unlist(mclapply(feats,filelines,mc.cores=cores))
		cat("pruning annotations\n")
		annos<-unlist(mclapply(annos,bedtoolsIntersect, b=targets, input="file", output="file", mc.cores=cores))
		feats<-unlist(mclapply(feats,bedtoolsIntersect, b=targets, input="file", output="file", mc.cores=cores))
	}
	
	cat("shuffling features\n")
	sfeats<-mclapply(1:numshuffles,function(z) unlist(lapply(feats,function(w) bed.shuffle(w,outname=paste0(tmpdir,"/",basename(removeext(w)),"_",z,"_shuffled.bed"),include=if(is.null(targetregions)==FALSE){targetregions} else{NULL}))) , mc.cores=cores , mc.preschedule=F)
	sfeatnames<-paste(featnames,"_shuffled",sep="")
	
	cat("counting features\n")
	numeachanno<-unlist(mclapply(annos,filelines,mc.cores=cores))
	numeachfeat<-unlist(mclapply(feats,filelines,mc.cores=cores))
	
	
	cat("calculating centers\n")
	if(usefeaturecenter==TRUE){
		feats<-unlist(mclapply(feats,bed.centers,mc.cores=cores))
		sfeats<-unlist(mclapply(sfeats,bed.centers,mc.cores=cores))
	}
	if(useannocenter==TRUE){
		annos<-unlist(mclapply(snnos,bed.centers,mc.cores=cores))
	}
	
	cat("finding nearby genes\n")
	if(is.null(genebed)==FALSE){
		feats<-unlist(mclapply(feats,bed.closest,bed2=genebed,mc.cores=cores))
		sfeats<-unlist(mclapply(sfeats,bed.closest,bed2=genebed,mc.cores=cores))
	}
	
	
	
	cat("finding overlaps\n")
	fxa<-as.data.frame(lapply(1:numannos, function(a){
		unlist(lapply(1:numfeats,function(f){
			as.numeric(system(paste("bedtools intersect -u -a",feats[f],"-b",annos[a],"| wc -l"),intern=TRUE))
			}))
		}))
	
	sfxa<-as.data.frame(mclapply(1:numannos, function(a){
		unlist(lapply(1:numfeats,function(f){
			mean(as.numeric(unlist(lapply(1:numshuffles, function(z){
				system(paste("bedtools intersect -u -a",sfeats[[z]][f],"-b",annos[a],"| wc -l"),intern=TRUE)
			}))))
		}))
	},mc.cores=cores,mc.preschedule=F))

	#bp<-matrix(numeachfeat,nrow=2,ncol=numannos)
	#bp<-sweep(bp,1,fxa,"-")

	#absolute #
	absolute<-t(
		matrix(
			c(
				unlist(as.data.frame(t(fxa))),
				rep(numeachfeat,each=numannos)-unlist(as.data.frame(t(fxa)))
			),
			nrow=numfeats*numannos , ncol=2)
	)
	fxa.pct<-100*fxa/numeachfeat
	sfxa.pct<-100*sfxa/numeachfeat
	fxa.ratio<-fxa.pct/sfxa.pct

s<-rep(0,numfeats*numannos)
s[(1:numfeats)*numannos+1]<-1
s<-s[1:(numfeats*numannos)]
barplot(unlist(as.data.frame(t(fxa.ratio))), col=rainbow(numannos),las=3,xaxt='n',space=s)
axis(1,(1:numfeats)*ceiling(numannos+1)-ceiling(numannos/2),labels=featnames,las=3)
legend("topleft",legend=annonames,fill=rainbow(numannos))


barplot(unlist(as.data.frame(t(fxa.pct))), col=rainbow(numannos),las=3,xaxt='n',space=s)




	#names=rep(annonames,numfeats) 
	barplot(absolute , names=rep(annonames,numfeats) , las=3)
	barplot(unlist(as.data.frame(t(fxa.pct))), col=rainbow(numannos),las=3,xaxt='n')
	axis(1,(1:numfeats)*numannos-numannos/2,labels=featnames)
	barplot(unlist(as.data.frame(t(sfxa.pct))), names=rep(annonames,numfeats) , las=3)
	barplot(unlist(as.data.frame(t(fxa.ratio))), names=rep(annonames,numfeats) , las=3)
	

	row.names(fxa.pct)=row.names(sfxa.pct)=row.names(fnxa.pct)=row.names(sfnxa.pct)=annonames
	colnames(fxa.pct)=featnames
	colnames(sfxa.pct)=sfeatnames
	colnames(fnxa.pct)=featnames
	colnames(sfnxa.pct)=sfeatnames
	a<-cbind(fxa.pct,sfxa.pct)
	b<-fxa.pct/sfxa.pct
	b[is.infinite(b)]<-NA
	pdf(file=paste("beddist_",suffix,".pdf",sep=""))
	par(mar=c(10,4,4,2))
	barplot(a,beside=T,las=2,ylab="% features overlapping with annotation",ylim=c(0,100),col=rainbow(numannos))
	legend("topright",legend=annonames,col=rainbow(numannos),lwd=3)
	
	bp<-barplot(b,beside=T,las=2,ylab="enrichment over shuffled",col=rainbow(numannos))
	bp<-barplot(b,beside=T,las=2,ylab="enrichment over shuffled",col=rainbow(numannos),ylim=c(0,bpylim))
	legend("topright",legend=annonames,col=rainbow(numannos),lwd=3)
	abline(h=1,col="grey30",lwd=1)
	text(x=bp,y=b,labels=as.numeric(b),xpd=TRUE,srt=90,adj=0)
	cat("check 1\n")
	
	
	if(is.null(scoremats) == FALSE & is.null(genebed)==FALSE){
		cat("pass1\n")
		matlist<-mclapply(scoremats,read.mat,mc.cores=cores)
		#if(do.call(identical,lapply(matlist,nrow))==FALSE){stop("matrices have different numbers of rows\n")}
		if(b73==TRUE){
			matgenes<-row.names(matlist[[1]])
		} else{matgenes<-unlist(lapply(strsplit(row.names(matlist[[1]]),"_"),"[",1))}
		
		fxa.genes<-lapply(1:numfeats, function(f) {
			lapply(1:numannos, function(a) {
				unique(readLines(pipe(paste("cut -f 4",fxa[a,f]))))
			})
		})
		sfxa.genes<-lapply(1:numfeats, function(f) {
			lapply(1:numannos, function(a) {
				unique(readLines(pipe(paste("cut -f 4",sfxa[a,f]))))
			})
		})
		fnxa.genes<-lapply(1:numfeats, function(f) {
			lapply(1:numannos, function(a) {
				unique(readLines(pipe(paste("cut -f 4",fnxa[a,f]))))
			})
		})
		
		sfnxa.genes<-lapply(1:numfeats, function(f) {
			lapply(1:numannos, function(a) {
				unique(readLines(pipe(paste("cut -f 4",sfnxa[a,f]))))
			})
		})
		
		
		all.genes<-unique(readLines(pipe(paste("cut -f 4",genebed))))
		all.featgenes<-lapply(1:numfeats,function(x) unique(readLines(pipe(paste("cut -f 4",feats[x])))))
		all.sfeatgenes<-lapply(1:numfeats,function(x) unique(readLines(pipe(paste("cut -f 4",sfeats[x])))))
		cat("checking if scoremats has values\n")
		for(s in 1:length(scoremats)){
			cat("finding gene scores for scoremats",s,"\n")
			fxa.genescores<-lapply(1:numfeats, function(f) {
				lapply(1:numannos, function(a) {
					if(fxa.counts[a,f]>0){
						rowMeans(as.matrix(matlist[[s]][which(matgenes %in% fxa.genes[[f]][[a]]),]),na.rm=TRUE)
					}
					else{
						NA
					}
				})
			})
			sfxa.genescores<-lapply(1:numfeats, function(f) {
				lapply(1:numannos, function(a) {
					if(sfxa.counts[a,f]>0){
						rowMeans(as.matrix(matlist[[s]][which(matgenes %in% sfxa.genes[[f]][[a]]),]),na.rm=TRUE)
					}
					else{
						NA
					}
				})
			})
			fnxa.genescores<-lapply(1:numfeats, function(f) {
				lapply(1:numannos, function(a) {
					if(fnxa.counts[a,f]>0){
						rowMeans(as.matrix(matlist[[s]][which(matgenes %in% fnxa.genes[[f]][[a]]),]),na.rm=TRUE)
					}
					else{
						NA
					}
				})
			})
			
			
			all.genescores<-rowMeans(matlist[[s]][which(matgenes %in% all.genes),],na.rm=TRUE)
			
			all.featgenescores<-lapply(1:numfeats,function(x) rowMeans(matlist[[s]][which(matgenes %in% all.featgenes[[x]]),],na.rm=TRUE))
			all.sfeatgenescores<-lapply(1:numfeats,function(x) rowMeans(as.matrix(matlist[[s]][which(matgenes %in% all.sfeatgenes[[x]]),]),na.rm=TRUE))
			
			
			#all<-list(1:numfeats, function(f) lapply(c(fxa.genescores,fnxa.genescores,sfxa.genescores), "[" ,  f ) )
			#all<-c(list(all.genescores),all.featgenescores,unlist(fxa.genescores,recursive=F),unlist(fnxa.genescores,recursive=F))
			
			all<-unlist(lapply(1:numfeats,function(f) unlist(lapply(1:numannos,function(a) list(fxa.genescores[[f]][[a]],fnxa.genescores[[f]][[a]],sfxa.genescores[[f]][[a]])),recursive=F)),recursive=F)
			#all<-unlist(lapply(1:numfeats,function(f) unlist(lapply(c(fxa.genescores,fnxa.genescores,sfxa.genescores), lapply , "[" ,  f ),recursive=F)),recursive=F)
			all<-c(list(all.genescores),all.featgenescores,all.sfeatgenescores,all)
			allnames<-basename(removeext(unlist(t(cbind(unlist(fxa),unlist(fnxa),unlist(sfxa))))))
			ats=c(1,1+(1:(numfeats*2)),unlist(lapply(1:(numfeats*numannos),function(x) x*3+1:3+x))+2*numfeats-2)
			par(mar=c(10,4,4,2))
			b<-boxplot(all,plot=F)
			boxplot(all,at=ats,col=c("white",rainbow(numfeats),rainbow(numfeats),rep(rainbow(numannos*numfeats),each=3)),names=c("all genes",featnames,featnames,allnames),las=3,cex.axis=0.7,outline=F,main=basename(removeext(scoremats[s])))
			#axis,1,at=ats,labels=
			points(ats,unlist(lapply(all,mean,na.rm=TRUE)))
			ats2<-b$stats[5,]
			text(ats,ats2,labels=paste(" ",lapply(all,length)),xpd=T,cex=0.7,srt=90,adj=0)
			
		}
		
		
	}
	dev.off()
}
