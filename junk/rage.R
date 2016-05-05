
# #######################################################################
#                R-based Analaysis of GEnomics data
#                              (rage)
#              written by Daniel Vera (dvera@fsu.edu)
# #######################################################################

cuffdiff.2.matrix <- function (cuffdiff.file , template.matrix.file , sample1name , sample2name , b73=FALSE){
	library(tools)
	matsuffix<-get.suffix(basename(template.matrix.file),"_")
	winsize=as.numeric(gsub("mat","",file_ext(template.matrix.file)))

	mat<-read.mat(template.matrix.file)
	cuf<-read.tsv(cuffdiff.file,header=T)
	matgenes<-unlist(lapply(strsplit(row.names(mat),";") , "[" , 3 ))
	if(b73){ matgenes<-remove.suffix(matgenes,"_T")}
	namematch<-match(matgenes,cuf[,3])
	 
	m1<-matrix(as.numeric(cuf[namematch,8]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(m1)<-row.names(mat)
	
	m2<-matrix(as.numeric(cuf[namematch,9]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(m2)<-row.names(mat)

	cuf$invq<-1/cuf[,13]
	gi2<-which(cuf[,9]-cuf[,8]>0)
	cuf$invq[gi2]<-cuf$invq[gi2]*-1
	invq<-matrix(as.numeric(cuf$invq[namematch]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(invq)<-row.names(mat)

	dif<-matrix(as.numeric(cuf[namematch,8]-cuf[namematch,9]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(dif)<-row.names(mat)

	invqname<-paste0(sample1name,"-",sample2name,"_inverseQval_",get.suffix(basename(template.matrix.file),"_"))
	m1name<-paste0(sample1name,"_FPKM_",matsuffix)
	m2name<-paste0(sample2name,"_FPKM_",matsuffix)
	difname<-paste0(sample1name,"Minus",sample2name,"_FPKM_",matsuffix)

	write.mat(m1,file=m1name)
	write.mat(m2,file=m2name)
	write.mat(dif,file=difname)
	write.mat(invq,file=invqname)
	
}
mat.tilegrid <- function (mat1 , mat2 , scoremats , mat1name = "x" , mat2name = "y" , LOG10=TRUE , abx = 0 , aby = 0 , legendnames = basename(removeext(scoremats)) , plotcolors=rainbow(length(scoremats)) , region = c(-1000,0) , tiles = 3 , ylims = c(-1,1) , breaks=NULL ){

	p=c(0,(1:(tiles-1))*1/tiles,1)
	s<-read.mat(mat1)
	r<-read.mat(mat2)
	sq<-quantile(s[,1],probs=p,na.rm=T)
	rq<-quantile(r[,1],probs=p,na.rm=T)
	cs<-cut(s[,1],sq,labels=F,include.lowest=T)
	cr<-cut(r[,1],rq,labels=F,include.lowest=T)
	
	sm2<-lapply(scoremats,read.mat)
	numscoremats<-length(scoremats)
	
	if(LOG10){
		ls=log10(s[,1]+1)
		lr=log10(r[,1]+1)
		xlabel=paste("log10 (",mat1name,")")
		ylabel=paste("log10 (",mat2name,")")
		lsq<-quantile(ls,probs=p,na.rm=T)
		lrq<-quantile(ls,probs=p,na.rm=T)
	} else{
		ls=s[,1]
		lr=r[,1]
		xlabel=mat1name
		ylabel=mat2name
		lsq<-sq
		lrq<-rq
	}

	plot(ls,lr,pch=16,xlab=xlabel,ylab=ylabel)
	abline(v=lsq[2:(tiles)],h=lrq[2:(tiles)],col="grey50")

	#X11()

	scatterdens(ls,lr,breaks=50)

	smoothScatter(ls,lr,xlab=xlabel,ylab=ylabel)
	abline(v=lsq[2:(tiles)],h=lrq[2:(tiles)],col="grey50")

	#X11()

	matcols<-ncol(sm2[[1]])
	xa<- ((0:4*(matcols/4)))
	ya<-c(ylims[1],(ylims[2]+ylims[1])/2,ylims[2])
	eg<-expand.grid(1:tiles,1:tiles)
	ml<-cbind(rbind(t(matrix(1:(tiles^2) , nrow=tiles)),(tiles^2+1):(tiles^2+tiles)),(tiles^2+tiles+1):((tiles+1)^2))
	ml<-ml[(tiles+1):1,]
	bottoms<-ml[tiles+1,]
	lefts<-ml[,1]
	layout(ml)
	par(oma=c(5,5,5,5),mar=c(0.2,0.2,0.2,0.2))

	for(i in 1:(tiles^2)){
		plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
		if(i %in% lefts ){
			axis(2,at=ya,labels=ya)
			mtext(which(i==rev(lefts)),side=2,line=2,cex=2)
		}
		if(i %in% bottoms){
			axis(1,at=xa,labels=xa)
			mtext(i,side=1,line=4,cex=2)
		}
		abline(v=abx,h=aby,col="grey70")
		for(j in 1:numscoremats){
			lines(1:matcols,colMeans(sm2[[j]][which(cs==eg[i,1] & cr==eg[i,2]),],na.rm=T),col=plotcolors[j])
		}
		text(xa[1],ya[3],paste0("n=",length(which(cs==eg[i,1] & cr==eg[i,2]))),pos=4)
	}

	for(i in 1:tiles){
		plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
		abline(v=abx,h=aby,col="grey70")
		if(i == 1){ axis(2,at=ya,labels=ya) }
		for(j in 1:numscoremats){
			lines(1:matcols,colMeans(sm2[[j]][which(cr==i),],na.rm=T),col=plotcolors[j])
		}
		text(xa[1],ya[3],paste0("n=",length(which(cr==i))),pos=4)
	}
	for(i in 1:tiles){
		plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
		abline(v=abx,h=aby,col="grey70")
		if(i == 1){ axis(1,at=xa,labels=xa) }
		for(j in 1:numscoremats){
			lines(1:matcols,colMeans(sm2[[j]][which(cs==i),],na.rm=T),col=plotcolors[j])
		}
		text(xa[1],ya[3],paste0("n=",length(which(cs==i))),pos=4)
	}
	plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
	abline(v=abx,h=aby,col="grey70")
	for(j in 1:numscoremats){
			lines(1:matcols,colMeans(sm2[[j]],na.rm=T),col=plotcolors[j])
	}
	text(xa[1],ya[3],paste0("n=",nrow(sm2[[1]])),pos=4)
	legend("topright",legend=legendnames , col=plotcolors , lwd=3 , bty="n")
	

}

coords.2.matrix <- function ( coordfile , genomefile , outdir , windowsize = 10000 , cores="max" ){
	
	library(parallel)
	library(gtools)
	if(cores=="max"){cores=detectCores()-1}
	
	cat("finding chromosomes\n")
	chroms <-  mixedsort (readLines ( pipe ( paste ( "cut -f 1" , coordfile , "| sort | uniq"))))
	chroms <- chroms[which(chroms != "chrM")]
	chroms <- chroms[which(chroms != "chrY")]
	numchroms <- length(chroms)
	chromsizes <- read.tsv ( genomefile )
	if(sum(chroms %in% chromsizes[,1]) != length(chroms)){stop("one or more chromosomes not found in genome file")}
	chromsizes<-chromsizes[match(chroms,chromsizes[,1]),]

	dir.create(outdir)
	outnames<-paste0(outdir,"/",chroms,"_win",windowsize,".mat")
		

	mclapply(1:numchroms,function(r) {
		cat("getting",chroms[r],"data\n")
		numchromwins<-floor(chromsizes[r,2]/windowsize)
		brks <- (0:numchromwins)*windowsize
		coords <- read.delim ( pipe ( paste0 ( "awk '($1==\"",chroms[r],"\" && $3==\"",chroms[r],"\")' ",coordfile," | cut -f 2,4" ) ) , stringsAsFactors=F , header = F )
		coords <- coords[which(coords[,1] <= brks[length(brks)] & coords[,2] <= brks[length(brks)] ),]
		a=coords[,1]
		b=coords[,2]
		
	
		cat("binning x-axis data\n")
		xbinind<-cut(a,brks,labels=F)
		
		cat("binning y-axis data\n")
		h<-mclapply(1:numchromwins,function(x){
			hist(b[which(xbinind==x)],breaks=brks,plot=F)
		},mc.cores=cores)
		
		cat("creating interaction matrix\n")
		binmat<-data.matrix(as.data.frame(lapply(h,"[[",2)))

		cat("saving binmat\n")
		write.tsv(binmat,file=outnames[r])
	} , mc.cores=cores , mc.preschedule=F )
}

matrix.2.tracks <- function ( matdir , genomefile , outname , windowsize = 10000 , maxdist = windowsize * 40 , cores="max" ){


	library(parallel)
	library(gtools)
	if(cores=="max"){cores=detectCores()-1}
	
	cat("finding chromosomes\n")
	matfiles<-mixedsort(files(paste0(matdir,"/*.mat")))
	chroms <-  unique(  ( remove.suffix( basename(matfiles),"_" ) ) )
	numchroms <- length(chroms)
	chromsizes <- read.tsv ( genomefile )
	if(sum(chroms %in% chromsizes[,1]) != length(chroms) ) {stop("one or more chromosomes not found in genome file")}
	chromsizes<-chromsizes[match(chroms,chromsizes[,1]),]
	numdistances=floor(maxdist/windowsize)
	d<-(1:numdistances)*windowsize
	distancenames<-formatC(d,width=nchar(maxdist),format="d",flag="0")

	outnames<-lapply(1:numchroms, function(o) paste0(matdir,"/",outname,"_",chroms[o],"_",distancenames,".bg"))
	foutnames<-paste0(matdir,"/",outname,"_",distancenames,".bg")
	
	for(b in 1:numchroms){
		

		
		mat<-read.tsv(matfiles[b])
		matcols<-ncol(mat)
		if(ceiling(matcols/2) != floor(matcols/2)){mat<-mat[1:(matcols-1),1:(matcols-1)]}
		matcols<-ncol(mat)
		l1<-lapply(1:matcols,function(x) 1:x)
		l2<-lapply(1:matcols,function(x) matcols+1-(x:1))
		rotmat<-matrix(NA,nrow=matcols,ncol=matcols)
		starts<-(matcols/2)-floor((0:(matcols-1)/2))
		
		cat("rotating",chroms[b],"matrix\n")
		fillers<-mclapply(1:matcols , function(i){
			v<-vector(length=i)
			for(j in 1:i){
				v[j]<-mat[l1[[i]][j],l2[[i]][j]]
			}
			return(v)
		} , mc.cores=cores , mc.preschedule=F)

		
		for(i in 1:matcols){
			rotmat[i,starts[i]:(starts[i]+i-1)]<-fillers[[i]]
		}
		
		rotmat[is.na(rotmat)]<-0

		oddstarts <- (0:(matcols-1))*windowsize + 1
		oddstops <- oddstarts + windowsize
		evenstarts <- oddstarts - windowsize/2
		evenstops <- evenstarts + windowsize
		md<-which(d>maxdist)[1]-1


		for(i in ((1:(numdistances/2))*2)-1){
			df<-data.frame(chroms[b],evenstarts,evenstops,rotmat[(matcols+1-i),])
			write.tsv(df,file=outnames[[b]][i])
		}

		for(i in ((1:(numdistances/2))*2)){
			df<-data.frame(chroms[b],oddstarts,oddstops,rotmat[(matcols+1-i),])
			write.tsv(df,file=outnames[[b]][i])
		}
	}

	for(b in 1:numdistances){
		infiles<-paste(unlist(lapply(outnames,"[",b)),collapse=" ")
		system(paste("cat",infiles,"| awk '($2 > 0)' >",foutnames[b]))
		system(paste("bedGraphToBigWig",foutnames[b],genomefile,gsub(".bg",".bw",foutnames[b])))
	}




}

replicationdomains<-function(bgfile , lspan = 0 , mergewithin = 100000 , removeY = T ){
	
	if(lspan != 0){ r<-bg.loess(r,lspan=lspan) }
	

	r<-read.tsv(bgfile)
	
	if(removeY==TRUE){
		ttr<-ttr[which(ttr$V1 != "chrY"),]
		r<-r[which(r$V1 != "chrY"),]
	}

	r$Vprod=c(1,r$V4[1:(rl-1)]*r$V4[2:rl])
	r$diff<-c(1,diff(r$V4))
	r$slope=c(1,diff(r$V4)/diff(r$V2))
	r$slopeprod=c(1,r$slope[1:(rl-1)]*r$slope[2:rl])
	spc<-which(r$slopeprod<0)
	pc<-which(r$prod<0)
	ttrb1<-unlist(lapply(1:length(pc), function(x) pc[x]-which(r[pc[x]:1,8]<0)[1]+1 ))
	ttrb2<-unlist(lapply(1:length(pc), function(x) pc[x]+which(r[pc[x]:pc[length(pc)],8]<0)[1]-1 ))
	
	chroms<-sort(unique(r$V1))
	numchroms<-length(chroms)
	chromlines<-lapply(1:numchroms,function(x) which(r[,1]==chroms[x]))

	
	t<-data.frame(pc,ttrb1,ttrb2)
	t[which(r[t[,2],2] > r[t[,1],2]),2]<-NA
	t[which(r[t[,3],2] < r[t[,1],2]),3]<-NA
	t$strand <- NA
	t[which(r[pc,4]<0),4]<-"+"
	t[which(r[pc,4]>0),4]<-"-"


	ttr<-data.frame("V1"=r[t[,1],1],"V2"=r[t[,2],2],"V3"=r[t[,3],3],"V4"=1:nrow(t),"V5"=1,"V6"=t[,4],"V7"=r[t[,1],2],"V8"=r[t[,1],3])
	ttr<-ttr[which(complete.cases(ttr)),]


	rd<-data.frame("V1"=r[t[1:(nrow(t)-1),1],1],"V2"=r[t[1:(nrow(t)-1),1],2],"V3"=r[t[2:(nrow(t)),1],2], "V4"=1:(nrow(t)-1) , "V5"=1 , "V6"=t[2:(nrow(t)),4] , "V7"=r[t[1:(nrow(t)-1),3],2],"V8"=r[t[2:(nrow(t)),2],3])
	otherchrom=r[t[2:(nrow(t)),1],1]
	rd<-rd[which(rd$V1==otherchrom),]
	rm(otherchrom)
	rd<-rd[which(complete.cases(rd)),]
}
bed.ends<-function( bedfile ){
	if(length(bedfile) > 1){stop("this function can only take 1 file. use lapply.")}
	outname<-paste0(basename(removeext(bedfile)),"_ends.bed")
	system(paste("awk '{print $1,$2,$2+1;print $1,$3,$3+1}' OFS='\t'",bedfile,"| sort --parallel=4 -T . -S 10G -k1,1 -k2,2n >",outname))
	return(outname)
}
fastq.clip.paired <- function( fastq1 , fastq2=NULL , adapter = "GATCGGAA" , quality.cutoff=20 , minlength = 10 , cores="max"){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	if(length(fastq1) < cores){cores=length(fastq1)}
	outnames<-paste0(basename(removeext(fastq1)),"_aclip.fastq")
	a<-mclapply(1:length(fastq1) , function(x) system(paste("cutadapt -O 1 -a",adapter,fastq1[x],">",outnames[x])) , mc.cores=cores , mc.preschedule=F)
	return(outnames)
}
hub.modify<- function(tracknames , hubloc , operation = "remove"){
	t<-grep("track ",a$V1)
	t2<-c(t[2:length(t)]-1,nrow(a))
	b<-strsplit(a$V1," ")
	bl<-unlist(lapply(b,length))
	d<-unlist(lapply(b,"[",1))
	e<-lapply(1:length(b), function(x) b[[x]][2:bl[x]] )
	e<-unlist(lapply(e,paste0,collapse=" "))
	f<-data.frame(V1=d,V2=e,stringsAsFactors=F)
	g<-lapply(1:length(t), function(x) data.frame(x=f[t[x]:t2[x],1] , y = f[t[x]:t2[x],2], stringsAsFactors=F))
	du<-unique(d)
	m<-lapply(1:length(g), function(x) match(g[[x]][,1],du))
	mat<-matrix(nrow=length(g),ncol=length(du))
	mat<-data.frame(mat,stringsAsFactors=F)
	colnames(mat) <- du
	for(i in 1:length(g)){
		mat[i,m[[i]]] = g[[i]][,2]
	}

	if(operation=="remove"){
		selected <- which(mat$track %in% tracknames | mat$parent %in% tracknames)
		if(length(selected) == 0 ){
			stop("no tracks found")
		} else{
			mat <- mat[-selected,]
			m <- m[-selected]
		}
	}

	if(operation=="duplicate"){
		selected <- which(mat$track %in% tracknames | mat$parent %in% tracknames)
		if(length(selected) == 0 ){
			stop("no tracks found")
		} else{
			mat <- mat[nrow(mat)+(1:length(selected)),] <- mat[selected,]
			m <- m[-selected]
		}
	}

	u<-unlist(
		lapply(
			1:nrow(mat), function(x) paste(colnames(mat)[m[[x]]],mat[x,m[[x]]])
		)
	)



	trackfile<-data.frame("V1"=u,stringsAsFactors=FALSE)
	write.tsv(trackfile,file="trackDb.txt")	

	# dont allow 'track' in track names
}
bgtocontigs<-function(bgfile,cores="max"){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	cat("loading bedgraph\n")
	bg<-read.tsv(bgfile)
	cat("splitting pos/neg\n")
	p<-bg
	n<-bg
	rm(bg)
	cat("converting scores\n")
	p$V4[which(p$V4<0)]<-0
	n$V4[which(n$V4>0)]<-0
	cat("finding gaps\n")
	d<-diff(p$V2)
	stops<-which(d !=1)
	starts=c(1,stops+1)
	stops=c(stops,nrow(p))
	dir.create("pos")
	dir.create("neg")
	dir.create("int")
	cat("saving data\n")
	a<-mclapply(1:length(starts),function(x){
		write.tsv(p[starts[x]:stops[x],4],file=paste("pos/",paste(p[starts[x],1:2],collapse="-"),".pos",sep=""))
		write.tsv(n[starts[x]:stops[x],4],file=paste("neg/",paste(n[starts[x],1:2],collapse="-"),".neg",sep=""))
		write.tsv(n[starts[x]:stops[x],1:3],file=paste("int/",paste(n[starts[x],1:2],collapse="-"),".int",sep=""))
	},mc.cores=cores)
	write.tsv(p,file=paste(basename(removeext(bgfile)),"_pos.bg",sep=""))
	write.tsv(n,file=paste(basename(removeext(bgfile)),"_neg.bg",sep=""))
}
bed.addcolor		<- function( bed, color, strand="." ){
	ext<-get.file.extensions(bed)
	curcol<-paste(col2rgb(color),collapse=",")
	curbed<-read.tsv(bed)
	outname<-paste0(removeext(bed),"_",color,".",ext)
	if(ncol(curbed)<3){ stop ("less than 3 columns in file")}
	if(ncol(curbed)==3){ curbed$V4 <- 1:nrow(curbed)}
	if(ncol(curbed)==4){ curbed$V5 <- 1}
	if(ncol(curbed)==5){ curbed$V6 <- strand}
	if(ncol(curbed)==6){ curbed$V7  <- curbed$V2}
	if(ncol(curbed)==7){ curbed$V8  <- curbed$V3}
	if(ncol(curbed)==8){ curbed$V9  <- curcol}
	write.tsv(curbed,file=outname)
}
bed.addstrand		<- function( bed, strand="+" ){
	ext<-get.file.extensions(bed)
	curbed<-read.tsv(bed)
	if(strand=="+"){ strandname=="plus" }
	if(strand=="-"){ strandname=="minus"}
	outname<-paste0(removeext(bed),"_",strandname,".",ext)
	if(ncol(curbed)<3){ stop ("less than 3 columns in file")}
	if(ncol(curbed)==3){ curbed$V4 <- 1:nrow(curbed)}
	if(ncol(curbed)==4){ curbed$V5 <- 1}
	if(ncol(curbed)==5){ curbed$V6 <- strand}
	write.tsv(curbed,file=outname)
}
stdout.tsv		<- function( cmd, header=FALSE, ... ){
	read.delim(pipe(cmd),header=header,stringsAsFactors=FALSE, ... )
}
matgrid<-function(mat){
	png(file="test.png",height=5000,width=1000)
	m<-matrix(1:3,nrow=3)
	mat[which(mat>5)]<-5
	mat[which(mat<1)]<-1
	rb<-colorRampPalette(c("black","red"))
	layout(m,heights=c(1,20,5),widths=1)
	par(oma=c(10,0,10,0))
	par(mar=c(10,15,5,15))
	image(matrix(1:100),col=rb(100),breaks=0:100,axes=F,bty="o")
	axis(side=1,at=c(0,1),labels=c(0,50),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
	mtext("sample1_0-200_scg3",side=3,cex=8,outer=T)
	par(mar=c(5,15,5,15),xpd=TRUE)
	image(t(mat),breaks=(1:50)/10,col=rb(49),axes=FALSE)
	lines(rep(-0.05,2),c(0,0.5),lwd=70,col="blue",lend=1)
	axis(side=1,at=c(0,0.5,1),labels=F,lwd=10,padj=1,line=1,tcl=-3)
	par(mar=c(20,15,0,15))
	plot(1:200,colMeans(mat,na.rm=T),type="l",lwd=10,xaxs="i",axes=FALSE)
	axis(side=1,at=c(1,100,200),labels=c(-1000,0,1000),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
	axis(side=2,cex.axis=8,lwd=10,tcl=-3,padj=-1)
	mtext("Distance from TSS (bp)",side=1,cex=5,outer=T,line=2)
	dev.off()
}
mat.plotmats		<- function( matlist , tracknames , colorlist="auto", plottypes = "l" , matnamelist="auto", ymaxs = "auto" , ymins = "auto" , cores = "max", centername="TSS" ){
	library(tools)
	library(parallel)
	if(cores=="max"){cores<-detectCores()-1}
	
	numsets<-length(matlist)
	nummats<-unlist(lapply(matlist,length))
	
	cat("loading matrices\n")
	mats<-lapply(1:numsets, function(s) lapply(1:nummats[s], function(m) read.mat( matlist[[s]][m] ) ) )
	
	if(matnamelist=="auto"){ matnamelist<-lapply(1:numsets, function(s) basename(removeext(matlist[[s]] ) ) ) }
	if(colorlist=="auto"){ colorlist<-lapply(1:numsets, function(s) rainbow(nummats[s]) ) }
	if(length(plottypes)==1){plottypes=rep(plottypes,numsets)}
	
	windowsize=as.numeric( gsub("mat","",file_ext(matlist[[1]][1]) ) )
	matcols<-ncol(mats[[1]][[1]])
	matrows<-nrow(mats[[1]][[1]])
	x<-((1:matcols)-(matcols/2))*windowsize
	xw<-c(min(x),((max(x)-min(x))*0.25+max(x)))
	
	
	cat("calculating y maxs\n")
	matmax<-lapply(1:numsets, function(s) lapply(1:nummats[s], function(m) apply(mats[[s]][[m]],1,max,na.rm=TRUE)  )  )
	yautomax<-lapply(1:numsets, function(s) unlist(mclapply(1:matrows, function(r) max(unlist(lapply(matmax[[s]],"[",r)), na.rm=TRUE) ,mc.cores=cores ) ) )
	
	cat("calculating y mins\n")
	matmin<-lapply(1:numsets, function(s) lapply(1:nummats[s], function(m) apply(mats[[s]][[m]],1,min,na.rm=TRUE)  )  )
	yautomin<-lapply(1:numsets, function(s) unlist(mclapply(1:matrows, function(r) min(unlist(lapply(matmin[[s]],"[",r)), na.rm=TRUE) ,mc.cores=cores ) ) )
	
	pdf(file="plotmatstest.pdf")
	par(mfrow=c(numsets,1),oma=c(4,0,4,0),mar=c(0,4,0,4))
	
	cat("drawing plots\n")
	for(r in 1:matrows){
		for(s in 1:numsets){
			plot(
				0,
				type="n",
				xaxt=if(s==numsets){'s'} else{'n'},
				xlim=xw,
				ylim=c(yautomin[[s]][r],yautomax[[s]][r]),
				ylab=tracknames[s],
				xlab=paste("Distance from",centername,"(bp)")
				#main=row.names(mats[[s]][[m]])[r]
			)
			grid(lty="solid",ny=NA)
			abline(v=max(x),lwd=4,col="black")
			for(m in 1:nummats[s]){
				lines(x,mats[[s]][[m]][r,],col=colorlist[[s]][m],type=plottypes[s])
			}
			legend("topright",legend=matnamelist[[s]],col=colorlist[[s]],lwd=3, cex=0.5, bty="n")
		}
	}
	dev.off()
}
isegtobed<-function(posfile,negfile=NULL,intervalfiles,prefix="z3b1sens",cutoffs=c(0,1,1.5,2,3),genome="b73v2", fix=FALSE){
	library(parallel)
	numcutoffs<-length(cutoffs)
	cutoffcols<-unlist(lapply(1:numcutoffs, function(x) 2+(x-1)*3 ))
	inames<-basename(removeext(intervalfiles))
	bothnames<-paste(prefix,"_BC",cutoffs,".bed",sep="")
	
	cat("loading files\n")
	intervals<-mclapply(intervalfiles,read.tsv,header=T,mc.cores=detectCores())
	pos<-read.tsv(posfile,header=T)
	posrows<-lapply(1:numcutoffs, function(x) which( is.na(pos[,cutoffcols[x]])==FALSE ) )
	posbgs<-lapply(1:numcutoffs, function(x) as.data.frame(matrix(NA,nrow=length(posrows[[x]]),ncol=6)))
	possubs<-lapply(1:numcutoffs, function(x) pos[posrows[[x]],])
	posmatches<-lapply(1:numcutoffs, function(x) match(possubs[[x]][,1],inames ) )
	posnames<-paste(prefix,"_pos_BC",cutoffs,".bed",sep="")
	allnames=posnames
	
	if(is.null(negfile)==FALSE){
		neg<-read.tsv(negfile,header=T)
		negrows<-lapply(1:numcutoffs, function(x) which( is.na(neg[,cutoffcols[x]])==FALSE ) )
		negbgs<-lapply(1:numcutoffs, function(x) as.data.frame(matrix(NA,nrow=length(negrows[[x]]),ncol=6)))
		negsubs<-lapply(1:numcutoffs, function(x) neg[negrows[[x]],])
		negmatches<-lapply(1:numcutoffs, function(x) match(negsubs[[x]][,1],inames ) )
		negnames<-paste(prefix,"_neg_BC",cutoffs,".bed",sep="")
		allnames<-c(posnames,negnames,bothnames)
	}
	
	
	for(i in 1:numcutoffs){
		cat("finding genomic locations of segments for BC",cutoffs[i],"\n")
		posbgs[[i]][,1]<-unlist(lapply(1:nrow(posbgs[[i]]), function(x) intervals[[posmatches[[i]][x]]][possubs[[i]][x,cutoffcols[i]],1] ))
		posbgs[[i]][,2]<-unlist(lapply(1:nrow(posbgs[[i]]), function(x) intervals[[posmatches[[i]][x]]][possubs[[i]][x,cutoffcols[i]],2] ))
		posbgs[[i]][,3]<-unlist(lapply(1:nrow(posbgs[[i]]), function(x) intervals[[posmatches[[i]][x]]][possubs[[i]][x,cutoffcols[i]+1]-1,3] ))
		#posbgs[[i]][,4]<-unlist(lapply(1:nrow(posbgs[[i]]), function(x) possubs[[i]][x,cutoffcols[i]+2] ))
		posbgs[[i]][,4]<-1:nrow(posbgs[[i]])
		posbgs[[i]][,5]<-1
		posbgs[[i]][,6]="+"
		
		
		badrows<-which(is.na(posbgs[[i]][,2]))
		if(length(badrows)>0){posbgs[[i]]<-posbgs[[i]][-badrows,]}
		
		write.tsv(posbgs[[i]],file=posnames[i])
		
		
		if(is.null(negfile)==FALSE){
			negbgs[[i]][,1]<-unlist(lapply(1:nrow(negbgs[[i]]), function(x) intervals[[negmatches[[i]][x]]][negsubs[[i]][x,cutoffcols[i]],1] ))
			negbgs[[i]][,2]<-unlist(lapply(1:nrow(negbgs[[i]]), function(x) intervals[[negmatches[[i]][x]]][negsubs[[i]][x,cutoffcols[i]],2] ))
			negbgs[[i]][,3]<-unlist(lapply(1:nrow(negbgs[[i]]), function(x) intervals[[negmatches[[i]][x]]][negsubs[[i]][x,cutoffcols[i]+1]-1,3] ))
			negbgs[[i]][,4]<-1:nrow(negbgs[[i]])
			negbgs[[i]][,5]<-1
			negbgs[[i]][,6]="-"
			badrows<-which(is.na(negbgs[[i]][,2]))
			if(length(badrows)>0){negbgs[[i]]<-negbgs[[i]][-badrows,]}
			write.tsv(negbgs[[i]],file=negnames[i])
			
			#both<-rbind(posbgs[[i]],negbgs[[i]],stringsAsFactors=FALSE)
			system(paste("cat",posnames[i],negnames[i],">",bothnames[i]))
			bed.sort(bothnames[i])
			#write.tsv(both,file=bothnames[i])
			
		}
	}
	
	allnames<-unlist(mclapply(allnames,bed.sort,mc.cores=detectCores()))
	allnames<-unlist(mclapply(allnames,bedToBigBed,genome=genome,mc.cores=detectCores()))
	return(allnames)
}
# ####################################
#      SYSTEM SETTINGS               #
# ####################################
hostname<-Sys.info()["nodename"]
user<-Sys.info()["user"]
if(hostname %in% c("genome","chrom1","chrom2","chrom3","chrom4","chrom5")){
    datapath<-"/data/"
}
if(grepl("spear",hostname)==TRUE){datapath<-"/lustre/maize/home/dlv04c/data/"}
if(grepl("submit",hostname)==TRUE){datapath<-"/lustre/maize/home/dlv04c/data/"}
if(grepl("medicine",hostname)==TRUE){datapath<-"/lustre/maize/home/dlv04c/data/"}
#use whole X11 terminal window
if(.Platform$GUI=="X11"){
	.adjustWidth <- function(...){
		options(width=Sys.getenv("COLUMNS"))
		TRUE
	}
}
if(.Platform$GUI=="X11"){
	.adjustWidthCallBack <- addTaskCallback(.adjustWidth)
}
options(scipen=100000000)
# ####################################
#      BASIC CONVENIENCE FUNCTIONS   #
# ####################################
`%ni%`			<- Negate(`%in%`)
rot			<- function( x){(1:x %% x) +1}
rotvec			<- function( vec){ vec[rot(length(vec))] }
rotrow			<- function( vec){ vec[rot(nrow(vec)),] }
rotcol			<- function( vec){ vec[,rot(ncol(vec))] }
rotlist			<- function( l ){lapply(c(2:length(l),1),function(x) l[[x]] ) }
lsl			<- function( ){system(command="ls -l")}
lsltr			<- function( ){system(command="ls -ltr")}
lslsr			<- function( ){system(command="ls -lSr")}
printname		<- function( df){
	print(deparse(substitute(df)))
}
remove.prefix		<- function( names,prefix){
	for(i in 1:length(names)){
		names[i]<-unlist(strsplit(names[i],prefix))[2]
	}
	names
}
get.prefix		<- function( names,separator){
	for(i in 1:length(names)){
		names[i]<-unlist(lapply(strsplit(names[i],separator),"[",1))
	}
}
get.suffix		<- function( names,separator){
	for(i in 1:length(names)){
		stringvec<-unlist(strsplit(names[i],separator))
		names[i]<-stringvec[length(stringvec)]
	}
	names
}
remove.suffix		<- function( names,suffix){
	for(i in 1:length(names)){
		names[i]<-unlist(strsplit(names[i],suffix))[1]
	}
	names
}
removeext		<- function( filenames ){
	filenames<-as.character(filenames)
	for(i in 1:length(filenames)){
		namevector<-unlist(strsplit(filenames[i],"\\."))
		filenames[i]<-paste(namevector[1:(length(namevector)-1)],collapse=".")
	}
	filenames
}
write.tsv		<- function( tsv, colnames=FALSE, rownames=FALSE, ... ){
	write.table(tsv,sep="\t",quote=FALSE,col.names=colnames,row.names=rownames, ... )
}
read.tsv			<- function( tsv, ... ){
	read.table(tsv,stringsAsFactors=FALSE,sep="\t", ... )
}
sc			<- function( cmd){
	system(command=paste(cmd))
}
findfiles		<- function( string, path="." ){
	return(readLines(pipe(paste("find ",path," -name '",string,"'",sep=""))))
}
files			<- function( x,...){
	parts<-(unlist(strsplit(x,"/")))
	#cat("parts",parts,"\n")
	if(length(parts)>1){
		pth<-paste(parts[1:length(parts)-1],collapse="/")
	}
	else{
		pth="."
	}
	#cat("pth",pth,"\n")
	nm<-parts[length(parts)]
	#cat("name",nm,"\n")
	if(length(parts)>1){
		list.files(pattern=glob2rx(nm),path=pth,full.names=TRUE, ... )
	}
	else{
		list.files(pattern=glob2rx(nm), ... )
	}
}
files2			<- function( cmd ){
	readLines(pipe(paste("ls -d",cmd)))
}
write.mat		<- function( mat, ... ){
	write.table(mat,sep="\t",quote=FALSE,col.names=FALSE, ... )
}
read.mat			<- function( mat, ... ){
	data.matrix(read.table(mat,stringsAsFactors=FALSE,sep="\t",row.names=1, ... ))
}
read.fmat		<- function( mat, ... ){
	as.matrix(read.table(mat,stringsAsFactors=FALSE,sep="\t",row.names=1, ... ))
}
downloadfile		<- function( url){
	if(Sys.info()["sysname"]=="Linux"){
		sc(paste("wget -q ",url,sep=""))
	}
	if(Sys.info()["sysname"]=="Darwin"){
		sc(paste("curl -silent -O ",url,sep=""))
	}
}
renamefiles		<- function( filelist, pattern="", replacement=""){
	if(pattern==""){stop("YOU MUST SPECIFY 'pattern'\n")}
	filenames<-basename(filelist)
	outnames<-gsub(pattern,replacement,filenames)
	nametab<-data.frame("before"=filenames,"after"=outnames,stringsAsFactors=FALSE)
	print(nametab)
	for(i in 1:length(filelist)){
		cat("renaming ",filenames[i]," to ",outnames[i],"\n",sep="")
		system(paste("mv",filelist[i],outnames[i]))
	}
}
filelines		<- function( filename ){
	as.numeric(readLines(pipe(paste("wc -l",filename,"| awk '{print $1}'"))))
}
getfullpaths		<- function( paths ){
	unlist(lapply(1:length(paths), function(x) system(paste("readlink -f",paths[x]),intern=TRUE) ))
}
uniquefilename		<- function( filename ){
	library(tools)
	if(grepl("\\.",basename(filename))==TRUE){
		ext<-paste(".",file_ext(filename),sep="")
	} else{ ext="" }
	suffixes<-c("",paste("_",1:100,sep=""))
	filenames<-paste(removeext(filename),suffixes,ext,sep="")
	filename<-filenames[which(file.exists(filenames)==FALSE)[1]]
	return(filename)
}
totalmem			<- function( ){
	return(as.numeric(readLines(pipe("free -b | awk '{print $2}' | head -n 2 | tail -n 1"))))
}
shead			<- function( filename, n=10 ){
	system(paste("head -n",n,filename))
}
stail			<- function( filename, n=10 ){
	system(paste("tail -n",n,filename))
}
# ####################################
#   FORMAT CONVERSION                #
# ####################################
sra.2.fastq		<- function( sras, paired=FALSE){
	library(parallel)
	arguments<-""
	if(paired==TRUE){arguments<-"--split-3"}
	doit<-mclapply(1:length(sras), function(x) system(paste("fastq-dump",arguments,sras[x])), mc.cores=detectCores() )
}
sam.2.bam		<- function( samfile, q=10 ){
	bamname<-paste(removeext(samfile),".bam",sep="")
	cat(bamname,": converting sam to bam\n")
	system(paste("samtools view -q",q,"-Sb ",samfile,">",bamname))
	return(bamname)
}
bam.2.bedpe		<- function( bamfile ){
	bamname<-basename(removeext(bamfile))
	outname<-paste(bamname,".bed",sep="")
	cat(bamname,": converting bam to bed\n")
	system(paste("bedtools bamtobed -bedpe -i",bamfile," | cut -f 1,2,6 | sort -V -T . -k1,1 -k2,2n -k3,3n > ",outname))
	return(outname)
}
bam.2.bed		<- function( bamfile ){
	bamname<-basename(removeext(bamfile))
	outname<-paste(bamname,".bed",sep="")
	cat(bamname,": converting bam to bed\n")
	system(paste("bedtools bamtobed -i",bamfile,">",outname))
	return(outname)
}
dz.2.wig			<- function( dzfile, cores="max", splitchromosomes=FALSE ){
	dzname<-basename(removeext(dzfile))
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	cat("finding chromosomes\n")
	chroms<-readLines(pipe(paste("cut -f 1",dzfile,"| uniq")))
	
	cat("splitting chromosomes\n")
	chromfiles<-unlist(mclapply(1:length(chroms), function(x){
		outname<-paste(dzname,"_",chroms[x],".wig",sep="")
		system(paste("grep -e $'",chroms[x],"\\t' ",dzfile," | cut -f 2,3 > ",outname,sep=""))
		system(paste("sed -i '1s/^/variableStep chrom=",chroms[x],"\\n/' ",outname,sep=""))
		return(outname)
	},mc.cores=cores,mc.preschedule=F))
	
	if(splitchromosomes==FALSE){
		outname<-paste(dzname,".wig",sep="")
		cat("concatenating chromosome scores\n")
		system(paste("cat",paste(chromfiles,collapse=" "),">",outname))
		system(paste("rm",paste(chromfiles,collapse=" ")))
	}
	else{ outname=chromfiles }
	return(outname)
}
scorelist.2.mat		<- function( scorefiles, genebed, namecol=1, delimiter="", namenumber=1, scorecol=4, nametype="symbol",matcol=200){
	
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
gff.2.bg			<- function( gff, extendbp=60, cores="max", genome="hg19"){
	library(gtools)
	cat("sorting gff\n")
	system(paste("sort -V -k1,1 -k4,4n",gff,"-o",gff))
	gffname<-basename(removeext(gff))
	outname<-paste(gffname,".bg",sep="")
	chroms<-readLines(pipe(paste("cut -f 1",getgenomefile(genome))))
	chroms<-chroms[mixedorder(chroms)]
	gff<-read.delim(pipe(paste("cut -f 1,4,5,6",gff)),stringsAsFactors=FALSE,header=FALSE)
	gff<-gff[which(gff$V1 %in% chroms),]
	if( identical(gff$V2,gff$V3)==TRUE ){
		cat(gffname,": start and stop coordinates identical, extending ends by",extendbp,"bp\n")
		gff$V3<-gff$V3+extendbp
	}
	cat(gffname,": saving to",outname,"\n")
	write.tsv(gff,file=outname)
	return(outname)
}
gff.2.bed			<- function( gff, extendbp=60, cores="max", genome="hg19", strand=TRUE , bedname=basename(removeext(gff)) ){
	library(gtools)
	cat("converting",gff,"\n")
	gffname<-basename(removeext(gff))
	outname<-paste(gffname,".bed",sep="")
	system(paste0("grep -v '#' ",gff," | awk 'BEGIN{i=0};{i+=1;print $1,$4,$5,\"",gffname,"_\"i,1,$7}' OFS='\t' > ",outname))




	# chroms<-readLines(pipe(paste("cut -f 1",getgenomefile(genome))))
	# chroms<-chroms[mixedorder(chroms)]
	# gff<-read.tsv(gff)
	# gff<-gff[order(gff$V1,gff$V4),]
	# gff<-gff[which(gff$V1 %in% chroms),]
	# if( identical(gff$V2,gff$V3)==TRUE ){
	# 	cat(gffname,": start and stop coordinates identical, extending ends by",extendbp,"bp\n")
	# 	gff$V3<-gff$V3+extendbp
	# }
	# gff<-data.frame(V1=gff$V1,V2=gff$V4,V3=gff$V5,V4=paste0(bedname,"_",1:nrow(gff)),V5=1,if(strand==T){V6=gff$V7},stringsAsFactors=F)
	# write.tsv(gff,file=outname)
	return(outname)
}
# ####################################
#   FASTQ FUNCTIONS                  #
# ####################################
fastq.fastqc		<- function( fastq ){
	system(paste("fastqc",fastq))
}
fastq.clip.fastx		<- function( fastq , adaptor="GATCGGAA", minlength=10  ){
	outname<-paste(basename(removeext(fastq)),"_clipped.fastq",sep="")
	system(paste("fastx_clipper -Q33 -v -n -k -l",minlength,"-a",adaptor,"-i",fastq,"-o",outname))
	return(outname)
}
fastq.clip		<- function( fastqs, paired=FALSE, filter=FALSE){
	r1names<-fastqs
	for(i in 1:length(r1names)){
	cat("trimming/filtering ",remove.suffix(r1names[i],"_R1.fastq"),"\n")
		r2name<-paste(remove.suffix(r1names[i],"_R1.fastq"),"_R2.fastq",sep="",collapse="")
		system(paste("/lustre/maize/home/dlv04c/software/trim_galore_zip/trim_galore --paired --length 10 --fastqc -q 25 ",r1names[i]," ",r2name,sep="",collapse=""))
		
	}
}
fastq.filter.fastx	<- function( fastq, qscore=20, percent=90, paired=FALSE){
	outname<-paste(basename(removeext(fastq)),"_filtered-q",qscore,"p",percent,".fastq",sep="")
	system(paste("fastq_quality_filter -Q33 -v -q",qscore,"-p",percent,"-i",fastq,"-o",outname))
	return(outname)
}
fastq.bt2		<- function( read1files, read2files=NULL,cores="max",genome="hg19",mode="end-to-end",input="-q", q=10, makebam = TRUE , makebed = TRUE , all=FALSE , dovetail=TRUE , discordant=FALSE , mixed=FALSE , unaligned=FALSE ){
	
	library(parallel)
	if(cores=="max"){cores<-detectCores()-1}
	
	indexfile<-getbt2index(genome)

	read1names<-basename(removeext(read1files))
	
	if(is.null(read2files)){
		paired=FALSE
		outnames<-paste(basename(removeext(read1files)),".sam",sep="")
	} else{
		paired=TRUE
		outnames<-paste(basename(removeext(read1files)),".sam",sep="")
	}
	
	for(i in 1:length(read1files)){
		cat(read1names[i],": aligning to genome\n")
		if(paired==TRUE){
			print(paste(
				"bowtie2",
				"-p",cores,
				if(dovetail){"--dovetail"},
				if(mixed==FALSE){"--no-mixed"},
				if(discordant==FALSE){"--no-discordant"},
				if(unaligned==FALSE){"--no-unal"},
				if(all){"-a"},
				input,
				"-x",indexfile,
				"-1",read1files[i],
				"-2",read2files[i],
				"-S",outnames[i]
			))
			system(paste(
				"bowtie2",
				"-p",cores,
				if(dovetail){"--dovetail"},
				if(mixed==FALSE){"--no-mixed"},
				if(discordant==FALSE){"--no-discordant"},
				if(unaligned==FALSE){"--no-unal"},
				if(all){"-a"},
				input,
				"-x",indexfile,
				"-1",read1files[i],
				"-2",read2files[i],
				"-S",outnames[i]
			))
		}
		else{
			system(paste(
				"bowtie2",
				"-p",cores,
				if(all){"-a"},
				input,
				"-x",indexfile,
				"-U",read1files[i],
				"-S",outnames[i]
			))
		}
	}
	
	if(cores > length(read1files)) { cores <- length(read1files) }
	if(makebam) { outnames<-unlist(mclapply(outnames,sam.2.bam,mc.cores=cores,q=q)) }
	if(makebed) { outnames<-unlist(mclapply(outnames,if(paired){bam.2.bedpe} else{bam.2.bed},mc.cores=cores))}

	return(outnames)
}


fastq.bt		<- function( read1files, read2files=NULL , cores="max" , genome="hg19" , input="-q", q=10, makebam = TRUE , makebed = TRUE , chunkmbs=64){
	
	library(parallel)
	if(cores=="max"){cores<-detectCores()-1}
	
	indexfile<-getbtindex(genome)

	read1names<-basename(removeext(read1files))
	
	if(is.null(read2files)){
		paired=FALSE
		outnames<-paste(basename(removeext(read1files)),".sam",sep="")
	} else{
		paired=TRUE
		outnames<-paste(basename(removeext(read1files)),".sam",sep="")
	}
	
	for(i in 1:length(read1files)){
		cat(read1names[i],": aligning to genome\n")
		if(paired==TRUE){
			system(paste(
				"bowtie",
				"-S",
				"--chunkmbs",chunkmbs,
				"-p",cores,
				input,
				indexfile,
				"-1",read1files[i],
				"-2",read2files[i],
				outnames[i]
			))
		}
		else{
			system(paste(
				"bowtie",
				"-S",
				"--chunkmbs",chunkmbs,
				"-p",cores,
				input,
				indexfile,
				read1files[i],
				outnames[i]
			))
		}
	}

	if(cores > length(read1files)) { cores <- length(read1files) }
	if(makebam) { outnames<-unlist(mclapply(outnames,sam.2.bam,mc.cores=cores,q=q)) }
	if(makebed) { outnames<-unlist(mclapply(outnames,if(paired){bam.2.bedpe} else{bam.2.bed},mc.cores=cores))}

	return(outnames)
}



fastq.bwasw		<- function( read1files, read2files=NULL , cores="max" , genome="hg19" , mode="end-to-end", input="-q", q=10, makebam = TRUE , makebed = TRUE ){
	
	library(parallel)
	if(cores=="max"){cores<-detectCores()-1}
	
	indexfile<-getbwaindex(genome)

	read1names<-basename(removeext(read1files))
	
	if(is.null(read2files)){
		paired=FALSE
		outnames<-paste(basename(removeext(read1files)),".sam",sep="")
	} else{
		paired=TRUE
		outnames<-paste(basename(removeext(read1files)),".sam",sep="")
	}
	
	for(i in 1:length(read1files)){
		cat(read1names[i],": aligning to genome\n")
		if(paired==TRUE){
			system(paste(
				"bwa bwasw",
				"-t",cores,
				indexfile,
				read1files[i],
				read2files[i],
				">",
				outnames[i]
			))
		}
		else{
			system(paste(
				"bwa bwasw",
				"-t",cores,
				indexfile,
				read1files[i],
				">",
				outnames[i]
			))
		}
	}
	
	if(cores > length(read1files)) { cores <- length(read1files) }
	if(makebam) { outnames<-unlist(mclapply(outnames,sam.2.bam,mc.cores=cores,q=q)) }
	if(makebed) { outnames<-unlist(mclapply(outnames,if(paired){bam.2.bedpe} else{bam.2.bed},mc.cores=cores))}

	return(outnames)
}




topcuff			<- function( read1file, read2file="", genome="hg19", librarytype="fr-firststrand", speed="--b2-very-fast", outdir=basename(removeext(read1file)) , cores="max"){
	library(parallel)
	if(cores=="max"){cores<-detectCores()-1}
	read1name<-basename(removeext(read1file))
	indexfile<-getbt2index(genome)
	tindexfile<-gettindex(genome)
	gtf<-getgtf(genome)
	cat(read1name,": aligning to genome\n")
	cat(paste("tophat -T -o ",outdir," --library-type ",librarytype," -p ",cores," --transcriptome-index=",tindexfile," ",indexfile," ",read1file," ",read2file,sep=""))
	cat("\n\n")
	system(paste("tophat -T -o ",outdir," --library-type ",librarytype," -p ",cores," --transcriptome-index=",tindexfile," ",indexfile," ",read1file," ",read2file,sep=""))
	bamfile<-(paste(outdir,"/accepted_hits.bam",sep=""))
	cat(paste("cufflinks","-G",gtf,"--library-type",librarytype,"-p",cores,"-o",outdir,bamfile))
	cat("\n\n")
	system(paste("cufflinks","-G",gtf,"--library-type",librarytype,"-p",cores,"-o",outdir,bamfile))
	return(paste(outdir,"/genes.fpkm_tracking",sep=""))
}
bam.sort			<- function( bamfile, memory=0.5*totalmem(), sortby="position" ){
	bamname<-basename(removeext(bamfile))
	#cat(bamname,": using",memory/1000000000,"Gb memory to sort bam by",sortby,"\n")
	if(sortby=="name"){extraargs="-n";suffix="_nsort"}
	if(sortby=="position"){extraargs="";suffix="_psort"}
	outname<-paste(bamname,suffix,sep="")
	#cat(paste("samtools sort",extraargs,bamfile,outname))
	
	system(paste("samtools sort",extraargs,bamfile,outname))
	return(outname)
}
bam.index		<- function( bamfile ){
	bamname<-basename(removeext(bamfile))
	cat(bamname,": indexing\n")
	system(paste("samtools index",bamfile))
}
# ####################################
#      BED FUNCTIONS                 #
# ####################################
cfbg.make		<- function( bedfile ){
	bedname<-basename(removeext(bedfile))
	outname<-paste(bedname,".cfbg",sep="")
	cat(bedname,": calculating fragment centers\n")
	system(paste("awk '{a=int(($2+$3)/2+0.5); $4=$3-$2; $2=a; $3=a+1;print}' OFS='\t' ",bedfile," | sort -V -T . -k1,1 -k2,2n > ",outname,sep=""))
	return(outname)
}
bed.wigaverages		<- function( beds, bigwigs, bednames=basename(removeext(beds)), wignames=basename(removeext(bigwigs)), targetregions="" ){
	numbeds<-length(beds)
	numwigs<-length(bigwigs)
	shufbeds<-unlist(lapply(beds,bed.shuffle,include=targetregions))
	shufbednames<-paste(bednames,"_shuffled",sep="")
	beds<-c(beds,shufbeds)
	bednames<-c(bednames,shufbednames)
	numbeds<-length(beds)
	bedlist<-lapply(beds,read.tsv)
	beds<-unlist(lapply(1:numbeds,function(x) {
		bed<-bedlist[[x]][,1:3]
		bed[,4]<-1:nrow(bed)
		outname<-paste(bednames[x],".tmp",sep="")
		write.tsv(bed,file=outname)
		outname
	}))

	scoresbybed<-lapply(1:numbeds,function(x){
		lapply(1:numwigs,function(y){
			outname<-paste(bednames[x],"_avg_",wignames[y],".tmp",sep="")
			system(paste("bigWigAverageOverBed",bigwigs[y],beds[x],outname))
			as.numeric(readLines(pipe(paste("cut -f 6",outname))))
		})
	})
	
	scoresbywig<-lapply(1:numwigs,function(x){
		lapply(lapply(scoresbybed,"[",x),unlist)
	})
	pdf(file="bedwigaverages.pdf")
	par(mar=c(10,4,4,2))
 	
	for(i in 1:numbeds){
		boxplot(scoresbybed[[i]],names=wignames,main=bednames[i],las=2,cex.names=0.5,outline=F)
	}
	for(i in 1:numwigs){
		boxplot(scoresbywig[[i]],names=bednames,main=wignames[i],las=2,cex.names=0.5,outline=F)
	}
	dev.off()
}
bed.shuffle		<- function( bedfile, include=NULL, exclude=NULL, samechrom=FALSE, chromfirst=FALSE, genome="hg19" , outname=NULL ){
	library(tools)
	bedname<-basename(removeext(bedfile))
	ext<-file_ext(bedfile)
	cat(bedname,": shuffling features\n")
	extraargs<-""
	if(samechrom==TRUE){extraargs<-"-chrom"}
	if(chromfirst==TRUE){extraargs<-paste(extraargs,"-chromFirst")}
	if(is.null(include) == FALSE){extraargs<-paste(extraargs,"-incl",include)}
	if(is.null(exclude) == FALSE){extraargs<-paste(extraargs,"-excl",exclude)}
	if(is.null(outname)){outname<-paste(bedname,"_shuffled.",ext,sep="")}
	print(paste("bedtools shuffle",extraargs,"-i",bedfile,"-g",getgenomefile(genome),">",outname))
	
	system(paste("bedtools shuffle",extraargs,"-i",bedfile,"-g",getgenomefile(genome),">",outname))
	return(outname)
}
bed.words		<- function( bedfiles, shufflewithin=NULL, genomefa="/data/hg19/igenome/Sequence/WholeGenomeFasta/genome.fa", numbases=50, sizerange=c(1,500), numfrags="all", reference="center", cores="max", strand=FALSE, lspan=0.07, slop=0 ){
	library(parallel)
	if(reference %in% c("center","end") == FALSE){stop("not a valid reference point, use 'center' or 'end'")}
	if(cores=="max"){cores=detectCores()-1}
	numbeds<-length(bedfiles)
	bednames<-basename(removeext(bedfiles))
		
	cat("determining dinucleotide words\n")
	nwords<-expand.grid(c("A","C","G","T"),c("A","C","G","T"))
	nnums<-expand.grid(1:4,1:4)
	dwords<-unique(paste(nwords[,1],nwords[,2],sep=""))
	numwords<-length(dwords)
	beds<-bedfiles
	beds<-unlist(mclapply(1:numbeds, function(i){
		
		cat(bednames[i],": processing\n")
		
		if(is.infinite(sizerange[2]) == FALSE & sizerange[1]>1){
			beds[i]<-bed.parselengths(beds[i],brks=sizerange)
		}

		if(numfrags != "all"){
			cat(bednames[i],": grabbing",numfrags,"random fragments\n")
			newname<-paste(bednames[i],"_",numfrags,"random.bed",sep="")
			system(paste("randomLines",beds[i],numfrags,newname))
			beds[i]<-newname
		}

		if(reference=="center"){
			cat(bednames[i],": calculating fragment centers\n")
			newname<-paste(bednames[i],"_centersforwords.bed",sep="")
			system(paste("awk '{a=int(($2+$3)/2);b=int(",numbases,"/2); $2=a-b-1-",slop,"; $3=a+b+1+",slop,";print}' OFS='\t' ",beds[i],"| sort -T . -k1,1 -k2,2n >",newname))
			beds[i]<-newname
		}
		
		return(beds[i])
	},mc.cores=detectCores()))
	
	
	
	basemats<-mclapply(1:numbeds, function(i){
		cat(bednames[i],": extracting sequences from bed\n")
		# FIX TO GET DIFFERENT GENOMES
		seqbedname<-paste(bednames[i],".seqbed",sep="")
		system(paste("bedtools getfasta -tab -fi",paste(genomefa,sep=""),"-bed",beds[i],"-fo",seqbedname))
		
		cat(bednames[i],": reading in sequences\n")
		basemat<-toupper(readLines(pipe(paste("cut -f 2",seqbedname))))
		basemat<-basemat[grep("N",basemat,invert=T)]
		numseqs<-length(basemat)
		
		cat(bednames[i],": creating matrix from sequences\n")
		basemat<-t(simplify2array(lapply(1:numseqs,function(x) substring(basemat[x],1:numbases,1:numbases))))
		if(strand==TRUE){
			cat("adjusting for strand\n")
			curbed<-read.tsv(beds[i])
			negrows<-which(curbed[,6] == "-")
			cat(length(negrows),"minus-stranded features found\n")
			basemat[negrows,1:numbases]<-basemat[negrows,numbases:1]
		}
		basemat
	},mc.cores=cores)
	
	factormats<-lapply(basemats,as.data.frame)
	nucfreqs<-mclapply(1:numbeds, function(i){
		freqs<-as.data.frame(lapply(1:ncol(factormats[[i]]),function(x) as.vector(table(factormats[[i]][,x]))))
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
	
	
	
	
	
	
	
	
	shufbeds<-unlist(mclapply(beds,bed.shuffle,include=shufflewithin))
	
	shufbednames<-paste(bednames,"_shuffled",sep="")
	shufbasemats<-mclapply(1:numbeds, function(i){
		cat(shufbednames[i],": extracting sequences from bed\n")
		# FIX TO GET DIFFERENT GENOMES
		seqbedname<-paste(shufbednames[i],".seqbed",sep="")
		system(paste("bedtools getfasta -tab -fi",paste(datapath,"hg19/igenome/Sequence/WholeGenomeFasta/genome.fa",sep=""),"-bed",shufbeds[i],"-fo",seqbedname))
		
		cat(shufbednames[i],": reading in sequences\n")
		basemat<-toupper(readLines(pipe(paste("cut -f 2",seqbedname))))
		numseqs<-length(basemat)
		

		cat(shufbednames[i],": creating matrix from sequences\n")
		basemat<-t(simplify2array(lapply(1:numseqs,function(x) substring(basemat[x],1:numbases,1:numbases))))
		if(strand==TRUE){
			cat("adjusting for strand\n")
			curbed<-read.tsv(shufbeds[i])
			negrows<-which(curbed[,6] == "-")
			cat(length(negrows),"minus-stranded features found\n")
			basemat[negrows,1:numbases]<-basemat[negrows,numbases:1]
		}
		basemat
	},mc.cores=cores)
	
	
	
	# ##################################
	#	SWITCH TO USING TABLE() ? MUST ADD ALL WORDS AT END TO MAKE SURE TABLE WORKS PROPERLY
	# #################################
	
	shuffreqmats<-list()
	for(j in 1:numbeds){

		cat(shufbednames[j],": finding dinucleotide frequencies\n")
		freqmat<-simplify2array(mclapply(1:(numbases-1), function(x) {
			curwords<-paste(shufbasemats[[j]][,x], shufbasemats[[j]][,x+1], sep="")
			unlist(lapply(1:numwords, function(y) length(which(curwords==dwords[y]))))
		},mc.cores=detectCores()))
		numseqs<-nrow(shufbasemats[[j]])
		#freqmat<-(freqmat/nrow(basemats[[j]]))*100
		freqmat<-rbind(freqmat,colSums(freqmat[c(1,4,13,16),]))
		freqmat<-rbind(freqmat,colSums(freqmat[c(6,7,10,11),]))
		row.names(freqmat)<-c(dwords,"TT/AA/TA/AT","CC/GG/CG/GC")
		shuffreqmats[[j]]<-freqmat
		#cat("saving dinucleotide frequency matrix\n")
		#write.mat(freqmat,file=paste(bednames[j],".dnfmat",sep=""))
		
	}
	logfreqmats<-freqmats
	logfreqmats<-lapply(1:numbeds,function(x) log2(freqmats[[x]]/shuffreqmats[[x]]))
	
	
	
	freqmats<-lapply(1:numbeds,function(x) (freqmats[[x]]/fragnums[x]))
	adjfreqmats<-lapply(1:numbeds,function(x) freqmats[[x]]/expectedwords[[x]])
	pdf(file="worddens.pdf")
	rb<-rainbow(numbeds)



	cat("plotting absolute nucleotide frequencies\n")
	for(j in 1:4){
		all<-unlist(lapply(1:numbeds,function(f) nucfreqs[[f]][j,]))
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(all),max(all)*1.2),main=row.names(nucfreqs[[1]]),xaxt='n',ylab="Frequency",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="black")
		}
		for(k in 1:numbeds){
			x<-1:(numbases)
			y<-unlist(nucfreqs[[k]][j,])
			if(lspan > 0){
				y<-loess(y~x,span=lspan)$fitted
			}
			lines(x,y,col=rb[k])
		}
		legend("topleft",legend=paste(bednames,"(",fragnums," fragments)"),col=rb,lwd=2)
		axis(side=1,at=seq(from=0,to=numbases,by=10))
	}


	for(j in 1:numbeds){
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(nucfreqs[[j]]),max(nucfreqs[[j]])),main=bednames[j],xaxt='n',ylab="Frequency",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="black")
		}
		lines(1:numbases,nucfreqs[[1]][1,],type="l",col="red")
		lines(1:numbases,nucfreqs[[1]][2,],type="l",col="green")
		lines(1:numbases,nucfreqs[[1]][3,],type="l",col="blue")
		lines(1:numbases,nucfreqs[[1]][4,],type="l",col="purple")
		legend("topright",legend=row.names(nucfreqs[[1]]),col=c("red","green","blue","purple"),lwd=3)
	}

	cat("plotting absolute dinucleotide frequencies\n")
	for(j in 1:18){
		all<-unlist(lapply(1:numbeds,function(f) freqmats[[f]][j,]))
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(all),max(all)*1.2),main=row.names(freqmat)[j],xaxt='n',ylab="Frequency (%)",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="black")
		}
		for(k in 1:numbeds){
			x<-1:(numbases-1)
			y<-freqmats[[k]][j,]
			if(lspan > 0){
				y<-loess(y~x,span=lspan)$fitted
			}
			#lines(x,y,type="s",col=rgb(0,0,0,75,maxColorValue=255)
			lines(x,y,col=rb[k])
		}
		legend("topleft",legend=paste(bednames,"(",fragnums," fragments)"),col=rb,lwd=2)
		axis(side=1,at=seq(from=0,to=numbases,by=10))
	}
	
	cat("plotting composite dinucleotide frequencies\n")
	
	for(j in 1:numbeds){
		all<-freqmats[[j]][c(17,18),]
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(all,na.rm=TRUE),max(all,na.rm=TRUE)*1.2),main=bednames[j],xaxt='n',ylab="Frequency (%)",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="black")
		}
		x<-1:(numbases-1)
		y1<-freqmats[[j]][17,]
		y2<-freqmats[[j]][18,]
		if(lspan > 0){
			y1<-loess(y1~x,span=lspan)$fitted
			y2<-loess(y2~x,span=lspan)$fitted
		}
		lines(x,y1,col="blue")
		lines(x,y2,col="red")
		legend("topleft",legend=c("AA/TT/AT/TA","CC/GG/CG/GC"),col=c("blue","red"),lwd=2)
		axis(side=1,at=seq(from=0,to=numbases,by=10))
	}
	
	




	cat("plotting adjusted dinucleotide frequencies\n")
	for(j in 1:18){
		all<-unlist(lapply(1:numbeds,function(f) adjfreqmats[[f]][j,]))
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(all),max(all)*1.2),main=row.names(adjfreqmats[[1]])[j],xaxt='n',ylab="log2 (observed / expected from nucleotide frequency)",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="grey50")
		}
		for(k in 1:numbeds){
			x<-1:(numbases-1)
			y<-unlist(adjfreqmats[[k]][j,])
			if(lspan > 0){
				y<-loess(y~x,span=lspan)$fitted
			}
			#lines(x,y,type="s",col=rgb(0,0,0,75,maxColorValue=255)
			lines(x,y,col=rb[k])
		}
		legend("topleft",legend=paste(bednames,"(",fragnums," fragments)"),col=rb,lwd=2)
		axis(side=1,at=seq(from=0,to=numbases,by=10))
	}
	
	cat("plotting adjusted composite frequencies\n")
	
	for(j in 1:numbeds){
		all<-adjfreqmats[[j]][c(17,18),]
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(all,na.rm=TRUE),max(all,na.rm=TRUE)*1.2),main=bednames[j],xaxt='n',ylab="log2 (observed / expected from nucleotide frequency)",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="grey50")
		}
		x<-1:(numbases-1)
		y1<-unlist(adjfreqmats[[j]][17,])
		y2<-unlist(adjfreqmats[[j]][18,])
		if(lspan > 0){
			y1<-loess(y1~x,span=lspan)$fitted
			y2<-loess(y2~x,span=lspan)$fitted
		}
		lines(x,y1,col="blue")
		lines(x,y2,col="red")
		legend("topleft",legend=c("AA/TT/AT/TA","CC/GG/CG/GC"),col=c("blue","red"),lwd=2)
		axis(side=1,at=seq(from=0,to=numbases,by=10))
	}
	cat("plotting normalized dinucleotide frequencies\n")
	
	for(j in 1:18){
		all<-unlist(lapply(1:numbeds,function(f) logfreqmats[[f]][j,]))
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(all,na.rm=TRUE),max(all,na.rm=TRUE)*1.2),main=row.names(freqmat)[j],xaxt='n',ylab="log2 ( observed / expected )",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="black")
		}
		for(k in 1:numbeds){
			x<-1:(numbases-1)
			y<-unlist(logfreqmats[[k]][j,])
			if(lspan > 0){
				y<-loess(y~x,span=lspan)$fitted
			}
			#lines(x,y,type="s",col=rgb(0,0,0,75,maxColorValue=255)
			lines(x,y,col=rb[k])
		}
		legend("topleft",legend=paste(bednames,"(",fragnums," fragments)"),col=rb,lwd=2)
		axis(side=1,at=seq(from=0,to=numbases,by=10))
	}
	
	
	cat("plotting composite normalized dinucleotide frequencies\n")
	
	for(j in 1:numbeds){
		all<-logfreqmats[[j]][c(17,18),]
		plot(0,type="n",xlim=c(0,numbases),ylim=c(min(all,na.rm=TRUE),max(all,na.rm=TRUE)*1.2),main=bednames[j],xaxt='n',ylab="log2 ( observed / expected )",xlab="nucleotide #")
		abline(v=seq(from=0,to=numbases,by=10),lwd=1,col="grey80")
		if(reference=="center"){
			abline(v=numbases/2,lwd=1,col="black")
		}
		x<-1:(numbases-1)
		y1<-logfreqmats[[j]][17,]
		y2<-logfreqmats[[j]][18,]
		if(lspan > 0){
			y1<-loess(y1~x,span=lspan)$fitted
			y2<-loess(y2~x,span=lspan)$fitted
		}
		#lines(x,y,typ
		#lines(x,y,type="s",col=rgb(0,0,0,75,maxColorValue=255)
		lines(x,y1,col="blue")
		lines(x,y2,col="red")
		legend("topleft",legend=c("AA/TT/AT/TA","CC/GG/CG/GC"),col=c("blue","red"),lwd=2)
		axis(side=1,at=seq(from=0,to=numbases,by=10))
	}
	
	
	dev.off()
}
bed.sample		<- function( bed , count ){
	totallines<-filelines(bed)
	if(count > totallines ){stop("sample quantity is greater than available lines")}
	outname<-paste(basename(removeext(bed)),"_random",count,".",get.file.extensions(bed),sep="")
	system(paste("randomLines",beds[i],numfrags,newname))
}
bed.spacing		<- function( beds, sizerange=c(1,500), distrange=c(1,1000), numfrags="all", legendnames=basename(removeext(beds)),plotcolors=rainbow(length(beds)),reference="center", cores="max", readsperblock=1000 ){
	library(parallel)
	if(reference %in% c("center","end") == FALSE){stop("not a valid reference point, use 'center' or 'end'")}
	if(cores=="max"){cores=detectCores()-1}
	numbeds<-length(beds)
	bednames<-basename(removeext(beds))
	
	
	beds<-unlist(mclapply(1:numbeds, function(i){
		
		if(is.infinite(sizerange[2]) == FALSE & sizerange[1]>1){
			beds[i]<-bed.parselengths(beds[i],brks=sizerange)
		}

		if(numfrags != "all"){
			cat(bednames[i],": grabbing",numfrags,"random fragments\n")
			newname<-paste(bednames[i],"_",numfrags,"random.bed",sep="")
			system(paste("randomLines",beds[i],numfrags,newname))
			beds[i]<-newname
		}

		if(reference=="center"){
			beds[i]<-bed.centers(beds[i])
		}
		beds[i]
	},mc.cores=detectCores()))
	
	
	bedcoords<-mclapply(1:numbeds, function(i){
		cat(bednames[i],"loading read coordinates\n")
		as.numeric(readLines(pipe(paste("cut -f 2",beds[i]))))
	},mc.cores=detectCores())
	
	dbedcoords<-mclapply(1:numbeds, function(i){
		cat(bednames[i],"loading duplicate-start read coordinates\n")
		as.numeric(readLines(pipe(paste("cut -f 1,2",beds[i],"| uniq -D | cut -f 2"))))
	},mc.cores=detectCores())
	
	distances<-lapply(1:numbeds, function(i){
		coords<-bedcoords[[i]]
		cat(bednames[i],"measuring distances for all reads\n")
		numblocks<-floor(length(coords)/1000)
		dh<-as.data.frame(mclapply(1:numblocks,function(x){
			d<-dist(coords[((((x-1)*1000):(x*1000))+1)])
			d<-na.omit(as.vector(d[d<=distrange[2] & d>=distrange[1] ]))
			h<-hist(d,breaks=distrange[1]:(distrange[2]+1),plot=FALSE)$counts/1000
			return(h)
		},mc.cores=cores))
		rowMeans(dh)
		
	})
	dmax<-quantile(unlist(distances),probs=0.98)
	
	ddistances<-lapply(1:numbeds, function(i){
		coords<-dbedcoords[[i]]
		cat(bednames[i],"measuring distances for all duplicate-start reads\n")
		numblocks<-floor(length(coords)/readsperblock)
		dh<-as.data.frame(mclapply(1:numblocks,function(x){
			d<-dist(coords[((((x-1)*readsperblock):(x*readsperblock))+1)])
			d<-na.omit(as.vector(d[d<=distrange[2] & d>=distrange[1] ]))
			h<-hist(d,breaks=distrange[1]:(distrange[2]+1),plot=FALSE)$counts/readsperblock
			return(h)
		},mc.cores=cores))
		rowMeans(dh)
		
	})
	ddmax<-quantile(unlist(ddistances),probs=0.98)

	cat("plotting all-read distogram\n")
	plot(0,type="n",xlim=c(distrange[1],distrange[2]),ylim=c(0,dmax),xlab="Distance (bp)",ylab="count", main="phasogram for all reads")
	abline(v=c((0:15)*10),lwd=1,col="grey80")
	for(i in 1:numbeds){
		lines(distrange[1]:distrange[2],distances[[i]],col=plotcolors[i],lwd=2)
	}
	legend("topright",legend=legendnames, col=plotcolors, lwd=3)
	
	X11()
	
	cat("plotting duplicate-start read distogram\n")
	plot(0,type="n",xlim=c(distrange[1],distrange[2]),ylim=c(0,ddmax),xlab="Distance (bp)",ylab="count", main="phasogram for duplicate-start reads")
	abline(v=c((0:15)*10),lwd=1,col="grey80")
	for(i in 1:numbeds){
		lines(distrange[1]:distrange[2],ddistances[[i]],col=plotcolors[i],lwd=2)
	}
	legend("topright",legend=legendnames, col=plotcolors, lwd=3)
}
bed.hist			<- function( beds , xlims=c(0,500), plotcolors=rainbow(length(beds)), legendnames=basename(removeext(beds)) ) {
	library(parallel)
	numbeds<-length(beds)
	
	cat("reading in beds\n")
	bedlist<-mclapply( 1:numbeds, function(x){
		as.numeric(readLines(pipe(paste("awk '{print $3-$2}'",beds[x]))))
	}, mc.cores=detectCores() )
	cat("calculating score distributions\n")
	densitylist<-mclapply( bedlist, density, from=xlims[1], to=xlims[2], mc.cores=detectCores() )
	numfrags<-unlist(lapply(bedlist, length)) 
	ymax<-max(unlist(lapply(1:numbeds, function(x) max(densitylist[[x]]$y ))))
	cat("plotting score distributions\n")
	plot(0,type="n",ylim=c(0,ymax), xlim=xlims, xlab="Fragment Size", ylab="Density")
	for(i in 1:numbeds){
		lines(densitylist[[i]],col=plotcolors[i],lwd=2)
	}
	legend("topright",legend=paste(legendnames,", n=",numfrags), col=plotcolors, lwd=3)
}
bed.structures <- function( bedfile , promoter5 = c(-1000,0) , promoter3 = c(0,1000) , genebody = c(200,-200) , bedname = basename(removeext(bedfile) ) ) {
	bed<-read.tsv(bedfile)
	bed<-bed[order(bed$V1,bed$V2),]
	if(ncol(bed) < 6){stop("no strand information found")}
	pos<-which(bed$V6=="+")
	neg<-which(bed$V6=="-")
	posbed<-bed[pos,]
	negbed<-bed[neg,]
	utr5<-rbind(
		data.frame(V1=posbed$V1,V2=posbed$V2,V3=posbed$V7,stringsAsFactors=F),
		data.frame(V1=negbed$V1,V2=negbed$V8,V3=negbed$V3,stringsAsFactors=F)
	)
	utr5<-utr5[which(utr5$V3-utr5$V2 > 0),]
	utr3<-rbind(
		data.frame(V1=posbed$V1,V2=posbed$V8,V3=posbed$V3,stringsAsFactors=F),
		data.frame(V1=negbed$V1,V2=negbed$V2,V3=negbed$V7,stringsAsFactors=F)
	)
	utr3<-utr3[which(utr3$V3-utr3$V2 > 0),]
	prom5<-rbind(
		data.frame(V1=posbed$V1,V2=posbed$V2+promoter5[1],V3=posbed$V2+promoter5[2],stringsAsFactors=F),
		data.frame(V1=negbed$V1,V2=negbed$V3-promoter5[2],V3=negbed$V3-promoter5[1],stringsAsFactors=F)
	)
	prom3<-rbind(
		data.frame(V1=posbed$V1,V2=posbed$V3+promoter3[1],V3=posbed$V3+promoter3[2],stringsAsFactors=F),
		data.frame(V1=negbed$V1,V2=negbed$V2-promoter3[2],V3=negbed$V2-promoter3[1],stringsAsFactors=F)
	)
	orf<-data.frame(V1=bed$V1,V2=bed$V7,V3=bed$V8,stringsAsFactors=F)
	gene<-data.frame(V1=bed$V1,V2=bed$V2+genebody[1],V3=bed$V3+genebody[2],stringsAsFactors=F)
	gene<-gene[which(gene$V3-gene$V2 > 0),]

	write.tsv(utr5,file=paste0(bedname,"_utr5.bed"))
	write.tsv(utr3,file=paste0(bedname,"_utr3.bed"))
	write.tsv(prom5,file=paste0(bedname,"_prom5.bed"))
	write.tsv(prom3,file=paste0(bedname,"_prom3.bed"))
	write.tsv(gene,file=paste0(bedname,"_genebody.bed"))
	write.tsv(orf,file=paste0(bedname,"_orf.bed"))


}
bed.dist			<- function( featurefiles,annotationfiles,suffix,numshuffles=100,bpylim=5,b73=FALSE,genebed=NULL, scoremats=NULL, targetregions=NULL, cores="max", featnames=basename(removeext(featurefiles)), annonames=basename(removeext(annotationfiles)),usefeaturecenter=FALSE,useannocenter=FALSE,plotcolors=rainbow(length(featurefiles))){
	
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
		targets<-bed.merge(targetregions)
		cat("counting total annotations\n")
		numtotalannos<-unlist(mclapply(annos,filelines,mc.cores=cores))
		numtotalfeats<-unlist(mclapply(feats,filelines,mc.cores=cores))
		cat("pruning annotations\n")
		annos<-unlist(mclapply(annos,bed.intersect, b=targets, input="file", output="file", mc.cores=cores))
		feats<-unlist(mclapply(feats,bed.intersect, b=targets, input="file", output="file", mc.cores=cores))
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
bed.recenter		<- function( bed,regionsize,center=FALSE,strand=TRUE,start=2,stop=3,input="object",mergebed=FALSE,buffer=0){
	if(input=="object"){
		bedname<-deparse(substitute(bed))
		windowbed<-bed
	}
	if(input=="file"){
		bedname<-basename(removeext(bed))
		windowbed<-read.tsv(bed)
	}
	if(strand==TRUE){
		negrows<-which(windowbed[,6]=="-")
		windowbed[negrows,2]<-windowbed[negrows,stop]
		windowbed[-negrows,2]<-windowbed[-negrows,start]
	}
	else{windowbed[,2]=windowbed[,start]}
	if(center==TRUE){
		windowbed[,2]=round((windowbed[,start]+windowbed[,stop])/2)
	}

	windowbed[,2]<-windowbed[,2]-regionsize/2
	windowbed[,3]<-windowbed[,2]+regionsize
	windowbed<-windowbed[which(windowbed[,2]> buffer),]
	winbedname<-paste(bedname,".winbed",sep="")
	write.tsv(windowbed,file=winbedname)
	return(winbedname)
}
bed.split		<- function( bed,regionsize,windowsize, output="file", input="object"){

	numwindows<-(regionsize/windowsize)
	
	if(input=="object"){
		bedname<-deparse(substitute(bed))
		curbed<-bed
	}
	if(input=="file"){
		bedname<-basename(removeext(bed))
		curbed<-read.tsv(bed)
	}
	
	
	bedrows<-nrow(curbed)
	bedcols<-ncol(curbed)
	
	#check parameters
	if(ceiling(numwindows)!=floor(numwindows)){stop("regionsize is not a multiple of windowsize")}
	if(ceiling(regionsize)!=floor(regionsize)){stop("regionsize must be an even number")}
	
	#make covbed
	flanks<- (  0:(numwindows-1)  ) * windowsize
	winstarts<-as.numeric( unlist( lapply( curbed[,2] , function(x) x + flanks ) ) )
	covbed<-data.frame("V1"=rep(curbed[,1],each=numwindows),"V2"=winstarts,"V3"=winstarts+windowsize,"V4"=rep(1:bedrows,each=numwindows))
	covbedname<-paste(bedname,".covbed",sep="")
	
	#return covbed
	write.tsv(covbed,file=covbedname)
	system(paste("sort -T . -S 10000000000b -k1,1 -k2,2n",covbedname,"-o",covbedname))
	return(covbedname)
}
bed.split.meta		<- function( bed, metasize, windowsize, flank, input="file", output="file", start=2, stop=3 ){
	library(parallel)
	if(input=="object"){
		bedname<-deparse(substitute(bed))
		curbed<-bed
	}
	if(input=="file"){
		bedname<-basename(removeext(bed))
		curbed<-read.tsv(bed)
	}

	bedrows<-nrow(curbed)
	bedcols<-ncol(curbed)

	numwindows<-metasize/windowsize
	numflankwindows<-flank/windowsize
	leftwinstarts<-0:(numflankwindows-1) * windowsize - flank
	leftwinends<-1:numflankwindows * windowsize - flank
	rightwinstarts<-0:(numflankwindows-1) * windowsize
	rightwinends<-1:numflankwindows * windowsize
	genesizes<-curbed[,stop]-curbed[,start]
	sizecos<-genesizes/metasize
	genewinstarts<-mclapply(1:bedrows, function(x) curbed[x,start] + 0:(numwindows-1) * windowsize * sizecos[x], mc.cores=detectCores())
	genewinends<-mclapply(1:bedrows, function(x) curbed[x,start] + 1:numwindows * windowsize * sizecos[x], mc.cores=detectCores())
	genewinstarts<-mclapply(genewinstarts,round,mc.cores=detectCores())
	genewinends<-mclapply(genewinends,round,mc.cores=detectCores())
	
	covbed<-data.frame(
		"V1"=rep(curbed[,1],each=numwindows+numflankwindows*2),
		"V2"=unlist(lapply(1:bedrows,function(x){ c(leftwinstarts+curbed[x,start],genewinstarts[[x]],rightwinstarts+curbed[x,stop]) } )),
		"V3"=unlist(lapply(1:bedrows,function(x){ c(leftwinends+curbed[x,start], genewinends[[x]],rightwinends+curbed[x,stop]) } )),
		"V4"=rep(1:bedrows,each=numwindows+numflankwindows*2)
	)
	covbedname<-paste(bedname,".covbed",sep="")
	#return covbed
	write.tsv(covbed,file=covbedname)
	system(paste("sort -T . -S 10000000000b -k1,1 -k2,2n",covbedname,"-o",covbedname))
	return(covbedname)
}
bed.sort			<- function( bedfile ){
	library(tools)
	ext<-file_ext(bedfile)
	system(paste("sort -T . -k1,1 -k2,2n",bedfile,"-o",bedfile))
	return(bedfile)
}
bed.parselengths		<- function( fragfile,brks=c(0,70,140,200),cores=1 ){
	
	#outputnames<-gsub("[[:punct:]]","",removeext(fragfiles))
	library(parallel)
	library(tools)
	if(cores=="max"){cores=detectCores()-1}
	if(length(fragfile) > 1){stop("bed.parselengths can only take 1 file")}
	numbreaks<-length(brks)
	brknames<-formatC(brks,width=nchar(brks[numbreaks]),format="d",flag="0")
	brknames[which(is.infinite(brks))]<-"Inf"
	print(brknames)
	ext<-file_ext(fragfile)
	fragname<-basename(removeext(fragfile))
	numfrags<-filelines(fragfile)
	scalar<-1000000/numfrags
	if(numbreaks==1){stop("need two or more breakpoints")}
	cat(fragname,": ",numfrags/1000000," million fragments\n",sep="")
	#parse fragments by length
	outnames<-unlist(mclapply(1:(numbreaks-1),function(i){
		cat(fragname,": ","processing fragments of length ",brks[i],"-",brks[i+1],"\n",sep="")
		subfragname<-paste0(fragname,"_",brknames[i],"-",brknames[i+1],".bed",sep="")
		if(brks[i+1]==Inf){
			print(paste("awk '$3-$2 >=",brks[i],"'",fragfile,">",subfragname))
			system(paste("awk '$3-$2 >=",brks[i],"'",fragfile,">",subfragname))
		} else{
			system(paste("awk '$3-$2 >=",brks[i],"&& $3-$2 <",brks[i+1],"'",fragfile,">",subfragname))
		}
		numsubfrags<-filelines(subfragname)
		cat(fragname,": ",numsubfrags,"/",numfrags," found (",100*numsubfrags/numfrags,"%)\n",sep="")
		return(subfragname)
	},mc.cores=cores))
	return(outnames)
}
bed.centers		<- function( bedfile ){
	if(length(bedfile) > 1){stop("bed.centers can only take 1 file")}
	#MAKE FILE OF FRAGMENT CENTERS
	outname<-paste(basename(removeext(bedfile)),"_centers.bed",sep="")
	cat("calculating fragment centers\n")
	system(paste("awk '{a=int(($2+$3)/2+0.5); $2=a; $3=a+1;print}' OFS='\t' ",bedfile," | sort -T . -k1,1 -k2,2n > ",outname,sep=""))
	return(outname)
}
bed.windowcov		<- function( bedfile, windowbed="", scalar="auto", windowsize=25, stepsize=windowsize, genome="hg19"){
	
	#check if only 1 file
	if(length(bedfile) > 1){stop("bed.windowcov can only take 1 file")}
	
	#get base name
	bedname<-basename(removeext(bedfile))
	outname<-paste(basename(removeext(bedfile)),"_win",windowsize,if(stepsize<windowsize){paste("_step",stepsize,sep="")},".bg",sep="")
	overlapsize<-(windowsize-stepsize)/2
	lsub=floor(overlapsize)
		
	#calculate scaling factor
	if(scalar=="auto"){
		cat(bedname,": calculating scalar\n")
		scalar=1000000/filelines(bedfile)
	}
	
	#make windows to calculate coverage over
	if(windowbed==""){
		cat(bedname,": making window file\n")
		wbedfile<-bedfile
		if(stepsize < windowsize){
			expand=round((windowsize-stepsize)/2)
			wbedfile<-bed.slop(wbedfile,expand,expand,genome=genome)
			wbedfile<-bed.merge(wbedfile,flank=expand)
		}
		windowbed=bed.makewindows(wbedfile,windowsize=windowsize,stepsize=stepsize)
	}
	#calculate coverage over windows
	cat(bedname,": calculating coverage\n")
	#system(paste("bedtools coverage -counts -a",bedfile,"-b",windowbed,"| awk '{print $1,$2,$3,$4*",scalar,"}' OFS='\\t' | sort -T . -k1,1 -k2,2n >",outname))
	
	system(paste(
		"bedtools coverage -counts -a",
		bedfile,
		"-b",
		windowbed,
		"| awk '{print $1,$2,$3,$4*",
		scalar,
		"}' OFS='\\t' |",
		if(stepsize<windowsize){paste("awk '{print $1,$2+",lsub,",$2+",lsub,"+",stepsize,",$4}' OFS='\t' |")},
		"sort -T . -k1,1 -k2,2n >",
		outname
	))
	
	
	cat(bedname,": making bigWig\n")
	outname<-bedGraphToBigWig(outname,genome=genome)
	return(outname)
}
bed.utrs			<- function( bedfile ){
}
bed.closest		<- function( bed1, bed2, strand=TRUE ){
	
	
	#ADD:
	# check if files exist
	# check if bed2 has a name column
	library(tools)
	ext<-file_ext(bed1)
	bedname<-basename(removeext(bed1))
	refname<-basename(removeext(bed2))
	curbed<-read.tsv(bed1)
	
	
	curbed$V4<-readLines(pipe(paste("cut -f 1,2,3",bed1,"| bedtools closest -t first -a /dev/stdin -b",bed2,"| cut -f 7")))
	if(strand==TRUE){
		curbed$V5<-1
		curbed$V6<-readLines(pipe(paste("cut -f 1,2,3",bed1,"| bedtools closest -t first -a /dev/stdin -b",bed2,"| cut -f 9")))
	}
	
	
	outname<-paste(bedname,"_closest-",refname,".",ext,sep="")
	write.tsv(curbed,file=outname)
	return(outname)
}
# ####################################
#      BEDTOOLS WRAPPER FUNCTIONS    #
# ####################################
bed.intersect		<- function( a, b, type="-u", fraction="", output="file", sorted=FALSE, input="file"){
	#TAKES (A) BED FILENAME AND KEEPS ONLY FEATURES THAT TOUCH (B), RETURNS RESULT OBJECT
	#if(input=="object"){
	#	bed1name<-deparse(substitute(a))
	#	bed2name<-deparse(substitute(b))
	#	write.tsv(paste(a,"*.bed",sep=""),file=bed1name)
	#	write.tsv(paste(b,"*.bed",sep=""),file=bed2name)
	#}
	
	
	extraargs<-""
	if(type != ""){extraargs<-type}
	if(fraction != ""){extraargs<-paste(extraargs,"-f",fraction)}
	if(sorted == TRUE){extraargs<-paste(extraargs,"-sorted")}
	
	
	if(input=="file"){
		bed1name<-basename(removeext(a))
		bed2name<-basename(removeext(b))
	}
	if(output=="object"){
		read.delim(pipe(paste("bedtools intersect ",extraargs," -a ",a," -b ",b,sep="")),header=FALSE,stringsAsFactors=FALSE)
	}
	if(output=="file"){
		filename<-paste(bed1name,"_x_",bed2name,".bed",sep="")
		system(paste("bedtools intersect ",extraargs," -a ",a," -b ",b," > ",filename,sep=""))
		return(filename)
	}
}
bed.merge		<- function( bedfile, flank=0, operation="none"){
    library(tools)
    ext<-file_ext(bedfile)
    if(flank==0){suffix=""}
    else{suffix=flank}
    if(operation=="none"){moreargs=""}
    else{
	suffix=c(suffix,operation)
	moreargs=c("-scores",operation)
    }
    
    outname<-paste(basename(removeext(bedfile)),"_merged",suffix,".",ext,sep="")
    
    system(paste("bedtools merge -d",flank,moreargs,"-i",bedfile,">",outname))
    return(outname)
}
bed.makewindows		<- function( bedfile, windowsize=25, stepsize=windowsize, mergebed=TRUE, mergeflank=500, outname="default", genome=FALSE){
	cat(bedfile,": making windows\n")
	if(outname=="default"){outname<-paste(basename(removeext(bedfile)),"_win",windowsize,".bed",sep="")}
	if(mergebed==TRUE & genome==FALSE){
		bedfile<-bed.merge(bedfile,flank=mergeflank)
	}
	if(windowsize>stepsize){
		extraargs<-"| awk 'x !~ $3; {x=$3}'"
	}
	else{
		extraargs<-""
	}
	if(genome==TRUE){
		bedfile<-getgenomefile(bedfile)
		inputarg="-g"
	} else{
		inputarg="-b"
	}
	print(paste("bedtools makewindows -w",windowsize,"-s",stepsize,inputarg,bedfile,extraargs,">",outname))
	
	system(paste("bedtools makewindows -w",windowsize,"-s",stepsize,inputarg,bedfile,extraargs,">",outname))
	return(outname)
}
bed.genomecov		<- function( bedfile, covmode="-bg" , genome="hg19", scalar="auto", makebigwig=TRUE, bam=FALSE ){
	if(length(bedfile) > 1){stop("bed.genomecov can only take 1 file")}
	if(scalar=="auto" & bam == FALSE){cat(bedfile,": counting fragments\n");scalar<-1000000/filelines(bedfile)}
	if(scalar=="auto" & bam == FALSE){cat(bedfile,": counting fragments\n");scalar<-1000000/bam.count(bedfile)}
	if(covmode=="-d"){pipes<-"| awk '{print $1,$2,$2+1,$3}' OFS='\t'"}
	if(covmode=="-dz"){pipes<-"| awk '{print $1,$2+1,$2+2,$3}' OFS='\t'"}
	if(covmode=="-bg"){pipes<-""}
	outname<-paste(basename(removeext(bedfile)),".bg",sep="")
	cat(bedfile,": calculating genome coverage\n")
	if(bam==TRUE){inputarg="-ibam"} else{inputarg="-i"}
	system(paste("bedtools genomecov",covmode,"-g",getgenomefile(genome),"-scale",scalar,inputarg,bedfile,pipes,">",outname))
	cat(bedfile,": genome coverage complete\n")
	if(makebigwig==TRUE){outname<-bedGraphToBigWig(outname,genome=genome)}
	return(outname)
}
bed.coverage		<- function( bedfile, covbedfile, numwindows, input="file", output="file"){
	if(length(bedfile) > 1){stop("bed.coverage can only take 1 file")}
	covlist<-read.delim(pipe(paste("bedtools coverage -counts -a ",bedfile," -b ",covbedfile," | sort -V -T . -k4,4n -k2,2n | cut -f 5",sep="")),header=FALSE,stringsAsFactors=FALSE)
	covlist<-t(matrix(covlist$V1,nrow=numwindows))
	covlist
}
bed.slop			<- function( bedfile, expandleft, expandright, genome="hg19" ){
	library(tools)
	bedname<-basename(removeext(bedfile))
	ext<-file_ext(bedfile)
	outname<-paste(basename(removeext(bedfile)),"_sl",expandleft,"_sr",expandright,".",ext,sep="")
	system(paste("bedtools slop","-i",bedfile,"-g",getgenomefile(genome),"-l",expandleft,"-r",expandright,">",outname))
	return(outname)
}
# ####################################
#      BEDGRAPH FUNCTIONS            #
# ####################################
bg.uniondiff		<- function( bglist1, bglist2, pattern="", replacement="", cores="max" ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	if(pattern=="" | replacement==""){stop("YOU MUST SPECIFY 'pattern' AND 'replacement' FOR OUTPUT FILE NAMES")}
	bglist1names<-basename(bglist1)
	bglist2names<-basename(bglist2)
	outbgnames<-gsub(pattern,replacement,bglist1names)
	nametab<-data.frame("BG1"=bglist1names,"BG2"=bglist2names,"OUTPUT"=outbgnames,stringsAsFactors=FALSE)
	print(nametab)
	if(length(which(bglist1names %in% outbgnames)==TRUE)>0){stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}
	cat("making unionbg\n")
	system(paste("bedtools unionbedg -i",paste(c(bglist1,bglist2),collapse=" "),"> tmpunionbg.tmp"))
	#cat(paste("bedtools unionbedg -i",paste(c(bglist1,bglist2),collapse=" "),"> tmpunionbg.tmp\n"))
	cat("calculating differences\n")
	mclapply(1:length(bglist1), function(x){
		#cat(paste("awk '{print $1,$2,$3,$",3+x,"-$",3+x+length(bglist1),"' OFS='\\t' tmpunionbg.tmp | cut -f 1,2,3,4 > ",outbgnames[x],"\n",sep=""))
		system(paste("awk '{print $1,$2,$3,$",3+x,"-$",3+x+length(bglist1),"}' OFS='\\t' tmpunionbg.tmp | cut -f 1,2,3,4 > ",outbgnames[x],sep=""))
	},mc.cores=cores)
}
bg.diff			<- function( bglist1, bglist2, operation="log2", pattern="", replacement="", cores="max" ){
	if(length(bglist1) != length(bglist2)){stop("number of files in bg lists don't match")}
	numbgs<-length(bglist1)
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	if(pattern=="" | replacement==""){stop("YOU MUST SPECIFY 'pattern' AND 'replacement' FOR OUTPUT FILE NAMES")}
	outbgnames<-gsub(pattern,replacement,basename((bglist1)))
	nametab<-data.frame("BG1"=bglist1,"BG2"=bglist2,"OUTPUT"=outbgnames,stringsAsFactors=FALSE)
	print(nametab)
	if(length(which(bglist1 %in% outbgnames)==TRUE)>0){stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}
	cat("loading files\n")
	bg1<-mclapply(bglist1,read.tsv,mc.cores=cores)
	bg2<-mclapply(bglist2,read.tsv,mc.cores=cores)
	bg1lines<-unlist(lapply(bg1,nrow))
	bg2lines<-unlist(lapply(bg2,nrow))
	if(bg1lines != bg2lines){stop("files do not match in rows")}
	cat("calculating differences\n")
	bg3<-mclapply(1:numbgs,function(x){
		if(operation=="difference"){
			bg1[[x]][,4]<-bg1[[x]][,4]-bg2[[x]][,4]
		}
		if(operation=="sum"){
			bg1[[x]][,4]<-bg1[[x]][,4]+bg2[[x]][,4]
		}
		if(operation=="log2"){
			bg1[[x]][,4]<-log2(bg1[[x]][,4]/bg2[[x]][,4])
		}
		if(operation=="mean" | operation=="average"){
			bg1[[x]][,4]<-(bg1[[x]][,4]+bg2[[x]][,4])/2
		}
		write.tsv(bg1[[x]],file=outbgnames[x])
	},mc.cores=cores)
	return(outbgnames)
}
bg.ops <- function(bglist , outputname , operation = "mean" , cores="max"){
	numbgs<-length(bglist)
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	bgs<-mclapply(bglist,read.tsv,mc.cores=cores)
	scores<-as.data.frame(lapply(bgs,"[",4))
	if(operation=="mean"){
		outscores<-rowMeans(scores,na.rm=T)
	}
	outbg<-data.frame(bgs[[1]])
	outbg$V4<-outscores
	write.tsv(outbg,file=outputname)

}
bg.qnorm			<- function( bgfiles, cores="max", comparable=TRUE ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	numfiles<-length(bgfiles)
	bgnames<-basename(removeext(bgfiles))
	outbgnames<-paste(bgnames,"_qnorm.bg",sep="")
	cat("loading files\n")
	bglist<-mclapply(bgfiles,read.tsv,mc.cores=cores)
	cat("determining if files are comparable\n")
	#if(unique(unlist(lapply(bglist,nrow))) > 1){comparable=FALSE}else{comparable=TRUE}
	if(comparable==FALSE){
		cat("making unionbg\n")
		allscores<-read.delim(pipe(paste("bedtools unionbedg -filler noscore -i ",paste(bgfiles,collapse=" ")," | grep -v 'noscore'",collapse=" ")),stringsAsFactors=FALSE,header=FALSE)
		
		cat("normalizing data\n")
		for(i in 2:numfiles){
			allscores[,3+i][order(allscores[,3+i])]<-allscores[,4][order(allscores[,4])]
		}
		
		cat("saving normalized data\n")
		for(i in 1:numfiles){
			write.tsv(cbind(allscores[,1:3],allscores[,3+i]), file=outbgnames[i])
		}
	}
	if(comparable==TRUE){
		cat("files are comparable\n")
		rankscores<-bglist[[1]][order(bglist[[1]][,4],na.last=F),4]
		bglist<-mclapply(1:numfiles, function(x){
			bglist[[x]][order(bglist[[x]][,4],na.last=F),4]<-rankscores
			write.tsv(bglist[[x]],file=outbgnames[x])
		},mc.cores=cores)
	}
}
bg.hist			<- function( bgfiles,X=TRUE,cores="max",xlims=NULL , legendnames=basename(removeext(bgfiles)),plotcolors=rainbow(length(bgfiles))){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	bgnames<-basename(removeext(bgfiles))
	numbgs<-length(bgfiles)
	densities<-mclapply(1:numbgs,function(x) density(as.numeric(readLines(pipe(paste("cut -f 4",bgfiles[x])))),na.rm=T),mc.cores=cores)
	ymax<-max(unlist(lapply(1:numbgs,function(x) densities[[x]]$y)))
	x<-unlist(lapply(1:numbgs,function(x) densities[[x]]$x))
	if(is.null(xlims)){xlims<-quantile(x,probs=c(0.05,0.95))}
	if(X==FALSE){pdf(file=paste(prefix,"_","densities.pdf",sep=""))}
	plot(0,type="n",xlim=xlims,ylim=c(0,ymax),ylab="Density",xlab="score")
	for(i in 1:numbgs){
		lines(densities[[i]],col=plotcolors[i],lwd=3)
	}
	legend("topleft",legend=legendnames, col=plotcolors, lwd=2)
	if(X==FALSE){dev.off()}
}
bg.dist			<- function( bgfile, bedfiles, cutoff = 1, piefeatures=c("intergenic","promoter5","5utr","CDS","3utr","promoter3"),cores="max" ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	numbeds<-length(bedfiles)
	bednames<-basename(removeext(bedfiles))
	
	cat("finding scores in annotations\n")
	annoscores<-mclapply(1:numbeds, function(a){
		as.numeric(readLines(pipe(paste("bedtools intersect -u -a",bgfile,"-b",bedfiles[a]," | cut -f 4"))))
	},mc.cores=cores)
	annoscores[[numbeds+1]]<-as.numeric(readLines(pipe(paste("cut -f 4",bgfile))))
	print(lapply(annoscores,length))
	
	bednames<-c(bednames,"all")
	numbeds<-numbeds+1
	
	#tabulate HS and HR bases
	cat("tabulating score summaries\n")
	scoredistmat<-matrix(ncol=numbeds,nrow=3)
	colnames(scoredistmat)<-bednames
	row.names(scoredistmat)<-c("neither","HS","HR")
	
	#counttotals<-unlist(mclapply(annoscores, length, mc.cores=cores) )
	#scoredistmat[2,]<-unlist(mclapply(1:numbeds,function(x) length(which(annoscores[[x]] >= cutoff)),mc.cores=cores ))
	#scoredistmat[3,]<-unlist(mclapply(1:numbeds,function(x) length(which(annoscores[[x]] <= -1*cutoff )),mc.cores=cores ))
	#scoredistmat[1,]<-counttotals-scoredistmat[2,]-scoredistmat[3,]
	#scoredistmat<sweep(scoredistmat,1, counttotals, "/")
	#scoredistmat<-scoredistmat[3:1,]
	#print(scoredistmat)
	
	
	pdf(file=paste(basename(removeext(bgfile)),"_cutoff",cutoff,"_annoscoredist.pdf",sep=""))
	
	cat("barplot\n")
	#barplot(scoredistmat,las=2,legend=rownames(scoredistmat),ylab="% genomic space")
	
	
	cat("density plot\n")
	rbc<-rainbow(numbeds)
	plot(0,type="n",ylim=c(0,1),xlim=c(-3,3),xlab="Light - Heavy Score",ylab="density")
	for(p in 1:numbeds){
		lines(density(annoscores[[p]],from=-3,to=3),col=rbc[p])
	}
	legend("topleft",legend=bednames, col=rbc, lwd=3)
	
	#cat("boxplot\n")
	#boxplot(annoscores,names=annos,las=2)
	
	cat("pie\n")
	pieannos<-which(bednames %in% piefeatures)
	pie(counttotals[pieannos],labels=paste(bednames[pieannos],  round(100*(counttotals[pieannos]))  /  sum(counttotals[pieannos])    ,"%"),clockwise=TRUE,main=paste("representation of each annotation out of the", sum(counttotals[pieannos])/1000000,"Mb"),col=rainbow(length(piefeatures)))
	
	dev.off()
	cat("finding max length of scores\n")
	maxlength<-max(unlist(mclapply(annoscores,length,mc.cores=detectCores())),na.rm=TRUE)
	cat("adjusting list lengths\n")
	annoscores<-mclapply(1:length(annoscores),function(x){ length(annoscores[[x]])=maxlength; annoscores[[x]]},mc.cores=detectCores())
	cat("making table of annoscores\n")
	annoscores<-do.call(cbind,annoscores)
	colnames(annoscores)<-bednames
	cat("saving table\n")
	write.table(annoscores,file="annoscores.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}
bg.window		<- function( bgfile, windowbed="", operation="sum" , mergebed=TRUE, genome="hg19", printzero=TRUE , windowsize=25, stepsize=windowsize, filler=0 ){
	bgname<-basename(removeext(bgfile))
	if(windowbed==""){
		cat(bgname,": windowing bedGraph regions\n")
		windowbed<-bed.makewindows(bgfile,windowsize=windowsize, stepsize=stepsize )
		outname<-paste(bgname,"_win",windowsize,".bg",sep="")
		prewin=FALSE
	}
	else{
		prewin=TRUE
		outname<-paste(bgname,"_win_",basename(removeext(windowbed)),".bg",sep="")
	}
	cat(bgname,": calculating",operation,"over window\n")
	
	windowedbg<-bg.map(bgfile,windowbed,operation=operation,filler=filler, outname=outname , printzero=printzero )
	if(stepsize < windowsize){
		overlapsize<-(windowsize-stepsize)/2
		lsub=floor(overlapsize)
		outname<-paste(removeext(outname),"_win",windowsize,"_step",stepsize,".bg",sep="")
		system(paste("awk '{print $1,$2+",lsub,",$2+",lsub,"+",stepsize,",$4}' OFS='\t' ",windowedbg," > ",outname,sep=""))
		windowedbg<-outname
	}
	cat(bgname,": making bigWig\n")
	outname<-bedGraphToBigWig(outname,genome=genome)
	return(windowedbg)
}
bg.loess			<- function( bgfile, lspan=0.05, cores="max" ){
	library(parallel)
	library(gtools)
	if(cores=="max"){cores=detectCores()-1}
	bgname<-basename(removeext(bgfile))
	cat(bgname,": loading file\n")
	curbg<-read.tsv(bgfile)
	curbg$V4[is.infinite(curbg$V4)]<-NA
	chroms<-unique(curbg$V1)
	chroms<-mixedsort(chroms)
	cat(bgname,": reference chromosome for span will be",chroms[1],"\n")
	numchroms<-length(chroms)
	pointsperchrom<-unlist(lapply(1:numchroms, function(x) nrow(curbg[which(curbg$V1==chroms[x]),])))
	outname<-paste(bgname,"_loess",gsub("\\.","-",lspan),".bg",sep="")
	cat(bgname,": smoothing chromosome data\n")
	lscores<-mclapply(1:numchroms, function(i){
		curchrom<-curbg[which(curbg$V1==chroms[i]),]
		if(i==1){clspan=lspan} else{clspan=lspan*(pointsperchrom[i]/pointsperchrom[1])}
		goodpoints<-which(complete.cases(curchrom[,4]))
		y<-curchrom[goodpoints,4]
		x<-curchrom[goodpoints,2]
		curchrom[goodpoints,4]<-loess(y~x,span=clspan)$fitted
		curchrom
	},mc.cores=cores,mc.preschedule=FALSE)
	curbg<-do.call(rbind,lscores)
	curbg<-curbg[order(curbg$V1,curbg$V2),]
	cat(bgname,": saving file\n")
	write.tsv(curbg,file=outname)
}
bg.plotchroms		<- function( bgfiles, outpdfname, overlay=TRUE, linewidth=2, ylabel="log2 ( early / late )" , ylims=c("auto","auto") , plotcolors=rainbow(length(bgfiles)) , legendnames = basename(removeext(bgfiles)) , cores="max", pdfdims=c(25,5)){
	library(parallel)
	library(gtools)
	if(cores=="max"){cores=detectCores()-1}
	numbgs<-length(bgfiles)
	bglist<-mclapply(bgfiles,read.tsv,mc.cores=cores)
	chroms<-mixedsort(unique(bglist[[1]][,1]))
	numchroms<-length(chroms)
	pdf(file=outpdfname, width=pdfdims[1], height=pdfdims[2])
	for(i in 1:numchroms){
		curchroms<-lapply(1:numbgs, function(z) bglist[[z]][which(bglist[[z]][,1]==chroms[i]),])
		chromsize<-max(curchroms[[1]][,3])
		allscores<-unlist(lapply(1:numbgs, function(z) curchroms[[z]][,4] ) )
		autoylims<-c(min(allscores,na.rm=TRUE),max(allscores,na.rm=TRUE))
		
		if(overlay==TRUE){
			plot(
				0,
				type="n",
				ylim=if(ylims[1]=="auto"){autoylims} else{ylims},
				xlim=c(0,chromsize),
				ylab=ylabel,
				xlab="chromosome coordinate (bp)",
				main=chroms[i],
				
			)
			abline(h=0,lwd=1,col="black")
			for(j in 1:numbgs){
				lines(curchroms[[j]][,2],curchroms[[j]][,4],lwd=linewidth, col=plotcolors[j])
			}
			legend("topright",legend=legendnames,col=plotcolors,lwd=4)
		}
		if(overlay==FALSE){
			par(mfrow=c(numbgs,1),oma=c(4,2,4,0),mar=c(0,4,0.75,10))
			for(j in 1:numbgs){
				plot(
					curchroms[[j]][,2],
					curchroms[[j]][,4],
					type="h",
					ylim=if(ylims[1]=="auto"){autoylims} else{ylims},
					xlim=c(0,chromsize),
					ylab="",
					xlab=if(j==numbgs){"chromosome coordinate (bp)"} else{""},
					xaxt=if(j==numbgs){'s'} else{'n'},
					col=plotcolors[j]
					
				)
				axis(side=1,labels=F)
				grid(lty="solid",ny=NA)
				abline(h=0,lwd=1,col="black")
				mtext(legendnames[j],side=4,las=1)
			}
			mtext(ylabel,side=2,outer=T,cex=1.5)
			mtext(chroms[i],side=3,outer=T,cex=1.5)
		}
	}
	dev.off()
}
bg.selectscores		<- function( bgfile, range=c(1,Inf) ){
	library(tools)
	bgname<-basename(removeext(bgfile))
	ext<-file_ext(bgfile)
	if(range[1]==-Inf){
		outname<-paste(bgname,"_lt",range[2],".",ext,sep="")
		system(paste("awk '$4 <=",range[2],"'",bgfile,">",outname))
	}
	if(range[2]==Inf){
		outname<-paste(bgname,"_gt",range[1],".",ext,sep="")
		system(paste("awk '$4 >=",range[1],"'",bgfile,">",outname))
	}
	if(range[1] != -Inf & range[2] != Inf){
		outname<-paste(bgname,"_gt",range[1],"_lt",range[2],".",ext,sep="")
		system(paste("awk '$4 >=",range[1],"'",bgfile,">",outname))
	}
	return(outname)
}
bg.threshold		<- function( bgfile, segmentnames=c("HS","HR"), segmergedist=75,cutoff=1, positive=TRUE, negative=FALSE){
	
	bgname<-basename(removeext(bgfile))
	
	cat("finding HS and HR regions\n")
	if(positive==TRUE){
		hsname<-paste(bgname,"_Min",cutoff,".bg",sep="")
		hs<-bg.selectscores(bgfile,range=c(cutoff,Inf))
		hs<-bed.merge(hs,flank=segmergedist)
		hs<-bg.map(bgfile,hs,operation="max",outname=hsname)
	}
	if(negative==TRUE){
		hrname<-paste(bgname,"_Max",cutoff,".bg",sep="")
		hr<-bg.selectscores(bgfile,range=c(-Inf,-1*cutoff))
		hr<-bed.merge(hr,flank=segmergedist)
		hr<-bg.map(bgfile,hr,operation="min",outname=paste(bgname,"_HR",cutoff,".bg",sep=""))
	}
	outnames<-c(if(positive==TRUE){hsname},if(negative==TRUE){hrname})
	return(outnames)
}
bg.averages		<- function( bedfile, covbedfile, numwindows, filler, input="file", output="file"){
	bedfile<-bed.sort(bedfile)
	#covbedfile<-bed.sort2(covbedfile)
	covlist<-read.delim(pipe(paste("bedtools map -c 4 -o mean -null \"",filler,"\" -a ",covbedfile," -b ",bedfile," | sort -T . -k4,4n -k2,2n | cut -f 5",sep="")),header=FALSE,stringsAsFactors=FALSE)
	covlist<-t(matrix(covlist$V1,nrow=numwindows))
	covlist
}
bg.map			<- function( bgfile, bedfile, operation="sum", filler=0, printzero = TRUE , outname="default" ){
    if(outname=="default"){outname<-paste(basename(removeext(bgfile)),"_map_",basename(removeext(bedfile)),"_",operation,".bg",sep="")}
    system(paste("bedtools map -c 4 -o",operation,"-null",filler,"-a",bedfile,"-b",bgfile,if(printzero==F){"| awk '{if ($4 != 0) print'}"},">",outname))
    return(outname)
}
bg.cors			<- function( bgfiles, outpdfname, bgnames=basename(removeext(bgfiles)), cores="max"  ){
	library(parallel)
	library(gplots)
	if(cores=="max"){cores=detectCores()-1}
	bgnames<-basename(removeext(bgfiles))
	numbgs<-length(bgfiles)
	
	colramp <- colorRampPalette(c("red","black","green"), space = "rgb")
	brks=seq(-1,1,by=2/100)
	hcolors<-colramp(100)
	
	#make list of matrices
	cat("reading in bedgraphs\n")
	bglist<-mclapply(bgfiles,read.tsv,mc.cores=cores)
	
	#make pairwise matrix correlations
	cat("calculating global correlations\n")
	pairs<-expand.grid(1:numbgs,1:numbgs)
	bgcors<-as.numeric(unlist(mclapply(1:nrow(pairs),function(x) cor(bglist[[pairs[x,1]]][,4] , bglist[[pairs[x,2]]][,4] , use="complete.obs" ),mc.cores=cores)))
	cormat<-matrix(bgcors,nrow=numbgs,ncol=numbgs)
	row.names(cormat)<-bgnames
	colnames(cormat)<-bgnames
	#for(i in 1:nrow(pairs)){
	#	cormat[pairs[i,1],pairs[i,2]]<-matcors[i]
	#}
	#globaltablename<-paste(matnames[1],"_globalcors.tsv",sep="")
	#cat("saving global correlation table to",globaltablename,"\n")
	#write.table(cormat,file=globaltablename,sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
	
	#plot pairwise global correlations as a heatmap
	cat("saving bg cross-correlation heatmap to",outpdfname,"\n")
	pdf(file=outpdfname)
	heatmap.2(cormat,Rowv=NA,Colv=NA,dendrogram="none",trace="none",col=hcolors,breaks=brks,main="Global cross-correlations")
	dev.off()
}
bg.transform		<- function( bgfiles,operation="unlog2",cores="max"){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	numbgs<-length(bgfiles)
	bgnames<-basename(removeext(bgfiles))
	outnames<-paste(bgnames,"_",operation,".bg",sep="")
	bgs<-mclapply(bgfiles,read.tsv,mc.cores=cores)
	bgs<-mclapply(1:numbgs,function(x){
		if(operation=="unlog2"){
			bgs[[x]][,4]=2^bgs[[x]][,4]
		}
		if(operation=="log2"){
			bgs[[x]][,4]=log2(bgs[[x]][,4])
		}
		if(operation=="sqrt"){
			bgs[[x]][,4]=sqrt(bgs[[x]][,4])
		}
		if(operation=="log10" | operation=="log"){
			bgs[[x]][,4]=log(bgs[[x]][,4])
		}
		if(operation=="unlog10"){
			bgs[[x]][,4]=10^bgs[[x]][,4]
		}
		write.tsv(bgs[[x]],file=outnames[x])
	},mc.cores=cores,mc.preschedule=FALSE)
	return(outnames)
}
# ####################################
#      VPLOT FUNCTIONS               #
# ####################################
vplot.make		<- function( fragmentfiles,features, regionsize=1000, windowsize=10, prunefeaturesto="", featurecenter=TRUE, strand=FALSE, start=2, stop=3, ylims=c(0,200), narrowpeak=FALSE){
	library(tools)
	library(parallel)
	numwindows<-regionsize/windowsize
	
	for(i in 1:length(features)){
		
		featname<-basename(removeext(features[i]))
		
		cat("\nstarting",features[i],"\n")
		if(narrowpeak==TRUE){   # ADD CODE TO CHECK IF REALLY NARROWPEAK (INSTEAD OF -1 IN COL 10)
			npbedname<-paste(basename(removeext(features[i])),".npbed",sep="")
			system(paste("awk '{$2=$2+$10; $3=$2+1;print}' OFS='\t'",features[i],">",npbedname))
			features[i]<-npbedname
		}
		
		numfeats<-filelines(features[i])
		cat(numfeats,"features\n")
		
		cat("making windowbed of features\n")
		features[i]<-bed.recenter(features[i],regionsize,center=featurecenter,strand=strand,start=start,stop=stop, input="file")
		
		if(prunefeaturesto != ""){
			cat("pruning features\n")
			features[i]<-bed.intersect(features[i],prunefeaturesto,input="file",output="file")
		}
		
		cat("reading in features\n")
		feats<-read.tsv(features[i])
		numfeats<-nrow(feats)
		cat(numfeats,"features\n")
		
		featnames<-feats$V4
		if(strand==TRUE){
			negrows<-which(feats[,6]=="-")
		}
		
		fragments<-fragmentfiles
		cat("making covbed\n")
		features[i]<-bed.split(features[i],regionsize,windowsize,input="file")
		
		for(j in 1:length(fragments)){
			fragname<-basename(removeext(fragments[j]))
			cat("\nstarting",fragname,"\n\t")
			
			
			
			if(grepl("tp://",fragments[j]) == TRUE ){
				cat("downloading file to current directory\n")
				system(paste("wget",fragments[j]))
				fragments[j]<-basename(fragments[j])
			}
			if(file_ext(fragments[j]) == "gz"){
				cat("extracting with gunzip to current directory\n")
				system(paste("gunzip -c",fragments[j],">",fragname ) )
				fragments[j]<-paste(basename(removeext(fragments[j])),sep="")
				fragname<-basename(removeext(fragments[j]))
				scoretype="bed"
			}
			if(file_ext(fragments[j]) %in% c("cbed","bed","cfbg","broadPeak","broadpeak","narrowPeak","narrowpeak") ==TRUE){
				scoretype="bed"
			}
			if(file_ext(fragments[j]) %in% c("bg","bedgraph","bedGraph") ==TRUE){
				scoretype="bedgraph"
			}
			if(file_ext(fragments[j]) %in% c("bb","bigbed","bigBed") == TRUE){
				cat("converting bigBed to bed\n")
				bigBedToBed(fragments[j])
				fragments[j]<-paste(fragname,".bed",sep="")
				scoretype="bed"
			}
			if(file_ext(fragments[j]) %in% c("bw","bigwig","bigWig") == TRUE){
				cat("converting bigWig to bedGraph\n")
				bigWigToBedGraph(fragments[j])
				fragments[j]<-paste(fragname,".bg",sep="")
				scoretype="bedgraph"
			}
			
			
			fragname<-basename(removeext(fragments[j]))
			
			
			
			
			vplotname<-paste(fragname,"_",featname,"_f",numfeats,".fmat",sep="")
			
			if(prunefragments==TRUE){
				cat("pruning fragments\n\t")
				fragments[j]<-bed.intersect(fragments[j],features[i],input="file",output="file")
			}
			
			cat("counting fragments...")
			numfrags<-filelines(fragments[j])
			cat(numfrags,"\n")
			



			cat("calculating fragment coverage around features\n\t")
			fmat<-read.delim(pipe(paste("bedtools map -c 4 -o collapse -null \"NA\" -a ",features[i]," -b ",fragments[j]," | sort -V -T . -k4,4 -k2,2n | cut -f 5",sep="")),header=FALSE,stringsAsFactors=FALSE)
			fmat<-t(matrix(fmat$V1,nrow=numwindows))
			
			if(strand==TRUE){
				cat("adjusting for strand\n")
				fmat[negrows,1:numwindows]<-fmat[negrows,numwindows:1]
			}
			
			#save fmat
			row.names(fmat)<-featnames
			write.mat(fmat,file=paste(basename(removeext(fragments[j])),"_",basename(featname),".fmat",sep=""))
			
			cat("processing fragments\n\t")
			fragsizes<-mclapply(1:numwindows,function(x) as.numeric(na.omit(as.numeric(unlist(strsplit(paste(fmat[,x],collapse=","),","))))),mc.cores=detectCores())
			
			cat("plotting fragment sizes in matrix\n\t")
			vmat<-matrix(0,ncol=numwindows,nrow=ylims[2]-ylims[1])
			for(h in 1:numwindows){
				vmat[,h]<-hist(fragsizes[[h]],breaks=ylims[1]:ylims[2],plot=F)$counts
			}
			
			#normalize counts
			if(rpm==TRUE){
				scalar<-1000000/numfrags
				vmat<-vmat*scalar
			}
			
			#vertically flip matrix
			vmat<-vmat[ylims[2]:1,]
			
			write.mat(,file=vplotname)
			vplot.draw(vplotname,png=TRUE)
			cat("vplot saved in",vplotname,"\n\n")
		}
	}
}
vplot.draw		<- function( vmat,hm.max="default",png=FALSE, plotcolors=c("black","green")){
	colramp <- colorRampPalette(plotcolors, space = "rgb")
	library(gplots)
	matname<-basename(removeext(vmat))
	m<-read.mat(vmat)
	if(hm.max=="default"){hm.max=as.numeric(round(quantile(m[m!=0],prob=0.99)))}
	filename<-paste(matname,"_max",round(hm.max),"_vplot.png",sep="")
	if(png==TRUE){
		png(file=filename,width=1000,height=1000)
	}
	heatmap.2(m,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=seq(0,hm.max,by=hm.max/100),labRow=NA,labCol=NA,col=colramp(100))
	if(png==TRUE){
		dev.off()
	}
}
vplot.rank		<- function( peakfiles, bgfiles, fragments, cellline="Gm12878", cores=16, np=TRUE, ... ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	cat("matching bedgraphs to peaks\n")
	bgfiles<-bgfiles[grep(cellline,bgfiles)]
	bgfiles<-bgfiles[grep("Rna",bgfiles,invert=T)]
	bgfiles<-bgfiles[grep("Tfbs",bgfiles)]
	peakfiles<-peakfiles[grep(cellline,peakfiles)]
	peakfiles<-peakfiles[grep("Rna",peakfiles,invert=T)]
	peakfiles<-peakfiles[grep("Tfbs",peakfiles)]
	bgnames<-basename(removeext(bgfiles))
	bgnames<-remove.prefix(bgnames,cellline)
	bgnames<-remove.suffix(bgnames,"Raw")
	bgnames<-remove.suffix(bgnames,"Sig")
	peaknames<-basename(removeext(peakfiles))
	peaknames<-remove.prefix(peaknames,cellline)
	peaknames<-remove.suffix(peaknames,"UniPk")
	peaknames<-remove.suffix(peaknames,"Pk")
	peaknames<-remove.suffix(peaknames,"Std")
	matches<-match(peaknames,bgnames)
	nobgs<-which(is.na(matches))
	matches<-na.omit(matches)
	cat(length(matches),"matches found\n")
	peakfiles<-peakfiles[-nobgs]
	bgfiles<-bgfiles[matches]
	numpeaks<-length(peakfiles)
	numbgs<-length(bgs)
	numfrags<-length(fragments)
	source("~/lus/mat.heatmap3.R")
	
	#for(p in 1:numpeaks){
	mclapply(1:numpeaks, function(p){
		
		dname<-mat.make( c(fragments, bgfiles[p]) , peakfiles[p], narrowpeak=np, fragmats=1:numfrags, cores=1, prunescores=T, featurecenter=T,prunefeaturesto="~/hg19/seqcap/seqcap_targets_merged.bed", maskbed="~/hg19/seqcap/seqcap_targets_merged.bed",regionsize=1000,closest="~/hg19/misc/scg3.bed" )
		hmbase<-basename(removeext(bgfiles[p]))
		hmname<-files(paste(dname,"/",hmbase,"*.mat10",sep=""))[1]
		hmnames<-files(paste(dname,"/*.mat10",sep=""))
		fmats<-files(paste(dname,"/*.fmat*",sep=""))
		mat.heatmap3( c(hmname,hmnames), fragmats=fmats,sorting=rep("mean,-50,50",2), plotcolors=c("deepskyblue black yellow"), cores=1 )
	},mc.cores=cores,mc.preschedule=F)
	#}
}
# ####################################
#      MATRIX FUNCTIONS              #
# ####################################
mat.make			<- function( scorefiles, features, closest=NULL, cores="max",meta=FALSE, metaflank=1000, maskbed=NULL,prunefeaturesto=NULL,rpm=TRUE,regionsize=2000,windowsize=10,start=2,stop=3,strand=TRUE,bgfiller=0,prunescores=FALSE,featurecenter=FALSE,suffix=NULL,scoremat=TRUE,fragmats=0, narrowpeak=FALSE){
        
        # TO DO
        # #######
        # add bed info to mat rownames
        # move tmp files to a tmp directory
        
        library(tools)
        #multicore setup
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	numwindows<-regionsize/windowsize
	numfeats<-length(features)
	
	#check parameters
	if(ceiling(numwindows)!=floor(numwindows)){stop("regionsize is not a multiple of windowsize")}
	if(ceiling(regionsize)!=floor(regionsize)){stop("regionsize must be an even number")}
	
	#process features
	for(i in 1:numfeats){
		
		featname<-basename(removeext(features[i]))
		
		#make directory for saving files
		dname<-paste(featname,"_mat",windowsize,sep="",collapse="")
		dname<-uniquefilename(dname)
		system(paste("mkdir",dname))
		
		if(is.null(closest) == FALSE){
			features[i]<-bed.closest(features[i],closest,strand=strand)
		}
		
		#recenter narrowpeaks
		if(file_ext(features[i]) %in% c("narrowPeak","narrowpeak","np") == TRUE & featurecenter==TRUE & narrowpeak==TRUE){
			cat("narrowPeak file detected\n")
			
			#check if peak location is in narrowPeak file and recenter features if OK
			headlines<-as.numeric(readLines(pipe(paste("head",features[i],"| awk '{print $10}'"))))
			if(length(which(headlines %in% (-1)>0))){
				cat("bad peak coordinate detected in narrowPeak file, treating as bed\n")
			} else{
				cat("peak coordinates detected in narrowPeak file, recentering bed features\n")
				outname<-paste(featname,"_recentered.np",sep="")
				system(paste("awk '{$2=$2+$10;$3=$2+$10+1;print}' OFS='\t' ",features[i]," | sort -T . -k1,1 -k2,2n > ",outname,sep=""))
				features[i]<-outname
			}
		}
		
		#create window bed for pruning
		cat("creating window bed\n")
		if(meta==FALSE){
			features[i]<-bed.recenter(features[i],regionsize,center=featurecenter,strand=strand,start=start,stop=stop,input="file" )
		}
		
		#prune features to specified file
		if(is.null(prunefeaturesto) == FALSE){
			cat("pruning features\n")
			features[i]<-bed.intersect(features[i],prunefeaturesto,input="file",output="file")
		}
		
		#copy processed features to output directory for later use
		system(paste("cp ",features[i]," ",dname,"/",sep=""))
		
		#read in bed and get dimensions
		cat("reading in bed and getting info\n")
		curbed<-read.tsv(features[i])
		
		if(meta==TRUE){
			curbed<-curbed[which(curbed[,2]>metaflank+1),]
			write.tsv(curbed,file=features[i])
		}
		
		bedcols<-ncol(curbed)
		bedrows<-nrow(curbed)
		
		# FEATURES SHOULD NOT BE MODIFIED AFTER THIS POINT #

		#find minus-stranded features
		if(bedcols < 6 & strand==TRUE){strand=FALSE;cat("WARNING: RUNNING WITH strand=FALSE BECAUSE NO STRAND COLUMN EXISTS\n")}
		if(strand==TRUE){negrows<-which(curbed[,6]=="-")}
		
		#make windows to calculate scores for matrix
		cat("making covbed\n")
		if(meta==FALSE){
			features[i]<-bed.split(features[i],regionsize,windowsize,input="file")
		} else{
			features[i]<-bed.split.meta(features[i],regionsize,windowsize, metaflank, input="file", start=start, stop=stop )
			numwindows<-(metaflank*2+regionsize)/windowsize
		}
		
		#make matrix of regions to (not) mask with NAs
		if(is.null(maskbed) == FALSE){
			cat("finding masking regions\n")
			maskmat<-bed.coverage(maskbed,features[i],numwindows)
			if(strand==TRUE){
				maskmat[negrows,1:numwindows]<-maskmat[negrows,numwindows:1]
			}
		}
		
		#assign row names to matrix
		cat("naming features\n")
		namebed<-curbed
		if(bedcols<4){
			namebed$V4<-1:bedrows
		}
		if(length(unique(namebed[,4])) != bedrows){
			namebed$V4<-paste(namebed$V4,1:bedrows,sep="-")
		}
		if(bedcols>12){
			colnames(namebed)[13]<-"symbol"
		} else{
			namebed$symbol<-namebed$V4
		}
		curbed[,4]<-paste(namebed$V4,namebed[,1],namebed$symbol,sep=";")
		rm(namebed)
		scores<-scorefiles
		

		
		#process each score file
		outs<-mclapply(1:length(scores), function(j) {
			
			scorename<-basename(removeext(scores[j]))
			
			#DOWNLOAD, EXTRACT, AND/OR FORMAT CONVERSION
			if(grepl("tp://",scores[j]) == TRUE ){
				cat(scorename,": downloading file to current directory\n")
				system(paste("wget",scores[j]))
				scores[j]<-basename(scores[j])
			}
			if(file_ext(scores[j]) == "gz"){
				cat(scorename,": extracting with gunzip to current directory\n")
				system(paste("gunzip -c",scores[j],">",scorename ) )
				scores[j]<-paste(basename(removeext(scores[j])),sep="")
				scorename<-basename(removeext(scores[j]))
			}
			if(file_ext(scores[j]) %in% c("wig","Wig") == TRUE){
				cat(scorename,": converting wig to bigWig\n")
				wigToBigWig(scores[j])
				scores[j]<-paste(scorename,".bw",sep="")
			}
			if(file_ext(scores[j]) %in% c("bb","bigbed","bigBed") == TRUE){
				cat(scorename,": converting bigBed to bed\n")
				bigBedToBed(scores[j])
				scores[j]<-paste(scorename,".bed",sep="")
			}
			if(file_ext(scores[j]) %in% c("bw","bigwig","bigWig") == TRUE){
				cat(scorename,": converting bigWig to bedGraph\n")
				bigWigToBedGraph(scores[j])
				scores[j]<-paste(scorename,".bg",sep="")
			}
			if(file_ext(scores[j]) %in% c("cbed","bed","cfbg","broadPeak","broadpeak","narrowPeak","narrowpeak") ==TRUE){
				scoretype="bed"
			}
			if(file_ext(scores[j]) %in% c("bg","bedgraph","bedGraph") ==TRUE){
				scoretype="bedgraph"
			}
			
			scorename<-basename(removeext(scores[j]))
			
			cat(scorename,": counting scores\n")
			numfrags<-filelines(scores[j])
			cat(scorename,":",numfrags,"fragments\n")
			
			#PRUNE READS/SCORES TO REGIONS AROUND FEATURES
			if(prunescores==TRUE){
				cat(scorename,": pruning scores to regions of interest\n")
				scores[j]<-bed.intersect(scores[j],features[i],input="file",output="file")
				cat(scorename,": counting pruned scores\n")
				prunedscorecount<-filelines(scores[j])
				cat(scorename,":",prunedscorecount,"scores after pruning\n")
			}
			
			if(scoremat == TRUE){
				#get coverage/average
				if(scoretype=="bed"){
					cat(scorename,": finding coverage of",scorename,"on",featname,"\n")
					curmat<-bed.coverage(scores[j],features[i],numwindows)
				}
				if(scoretype=="bedgraph"){
					cat(scorename,": mapping scores in",scorename,"to",featname,"\n")
					curmat<-bg.averages(scores[j],features[i],numwindows,bgfiller)
				}
				
				#add row names to matrix
				row.names(curmat)<-curbed[,4]
				
				#flip rows of negative-stranded features
				if(strand==TRUE){
					curmat[negrows,1:numwindows]<-curmat[negrows,numwindows:1]
				}
				
				#normalize by total fragments
				if(rpm==TRUE & scoretype=="bed"){
					cat(scorename,": normalizing data\n")
					scalar<-1000000/numfrags
					curmat<-curmat*scalar
				}
				
				#mask uncovered regions with NA and remove poorly-covered regions
				if(is.null(maskbed) == FALSE){
					cat(scorename,": masking data\n")
					curmat[maskmat==0]<-NA
				}
				
				#save matrix
				outfilename<-paste(scorename,"_",featname,suffix,".mat",windowsize,sep="")
				write.mat(curmat,file=paste(dname,"/",outfilename,sep=""))
				cat(scorename,": matrix saved to",paste(dname,"/",outfilename,sep=""),"\n")
				return(paste(dname,"/",outfilename,sep=""))
			}
			return(outs)
			if(j %in% fragmats & scoretype == "bed"){
				cat("boing\n")
				cfbg<-cfbg.make(scores[j])
				cat(scorename,": calculating fragment coverage around features for fragmat\n")
				
				fmat<-readLines(pipe(paste("bedtools map -c 4 -o collapse -null \"NA\" -a ",features[i]," -b ",cfbg," | sort -V -T . -k4,4 -k2,2n | cut -f 5",sep="")))
				fmat<-t(matrix(fmat,nrow=numwindows))
				
				if(strand==TRUE){
					cat("adjusting for strand\n")
					fmat[negrows,1:numwindows]<-fmat[negrows,numwindows:1]
				}
				print(dim(fmat))
				#save fmat
				row.names(fmat)<-curbed[,4]
				#headername<-paste(removeext(fmatname),".header",sep="")
				#tmpname<-paste(scorename,"_",featname,suffix,".tmp",sep="")
				#fmatname<-paste(scorename,"_",featname,suffix,".fmat",windowsize,sep="")
				fmatname<-paste(scorename,"_",featname,suffix,".fmat",windowsize,sep="")
				cat("saving matrix to",fmatname,"\n")
				write.mat(fmat,file=paste(dname,"/",fmatname,sep=""))
				#cat("adding header to",fmatname,"\n")
				#system(paste("echo '#",numfrags,"' > ",headername," && cat ",headername," ",tmpname," > ",fmatname,sep=""))
			}
		#})
		},mc.cores=cores, mc.preschedule=FALSE)
	}
	
}
mat.diff			<- function( matlist1,matlist2,operation="subtract",pattern="",replacement=""){
	if(pattern=="" | replacement==""){stop("YOU MUST SPECIFY 'pattern' AND 'replacement' FOR OUTPUT FILE NAMES")}
	matlist1names<-basename(matlist1)
	matlist2names<-basename(matlist2)
	outmatnames<-gsub(pattern,replacement,matlist1names)
	nametab<-data.frame("MATRIX1"=matlist1names,"MATRIX2"=matlist2names,"OUTPUT"=outmatnames,stringsAsFactors=FALSE)
	print(nametab)
	if(length(which(matlist1names %in% outmatnames)==TRUE)>0){stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}
	for(i in 1:length(matlist1)){
		cat("processing ",outmatnames[i],"\n",sep="")
		mat1<-read.mat(matlist1[i])
		mat2<-read.mat(matlist2[i])
		if(operation=="subtract"){
			mat3<-mat1-mat2
		}
		if(operation=="log2"){
			mat3<-log2(mat1/mat2)
			mat3[is.infinite(mat3)]<-NA
		}
		write.mat(mat3,file=outmatnames[i])
	}
}
mat.cors			<- function( mats, numgroups=3,hcluster=FALSE,plotcolors=c("red","green","black"),sorton=2,legendnames=basename(removeext(mats))){
	library(parallel)
	library(gplots)
	colramp <- colorRampPalette(plotcolors, space = "rgb")
	brks=seq(-1,1,by=2/100)
	hcolors<-colramp(100)
	matnames<-basename(removeext(mats))
	nummats<-length(mats)
	
	#make list of matrices
	cat("reading in matrices\n")
	matlist<-mclapply(1:nummats,function(x){  read.mat(mats[x]) },mc.cores=detectCores() )
	numcols<-ncol(matlist[[1]])
	#remove rows with no variance or no scores
	cat("finding genes with no variance\n")
	badvarrows<-unique(unlist(lapply(1:nummats, function(x){ which(apply(matlist[[x]],1,sd, na.rm=TRUE)==0) } )))
	
	cat("finding genes with no scores\n")
	badrows<-unique(unlist(lapply(1:nummats, function(x){ which(rowSums(is.na(matlist[[x]]))==numcols) })))
	
	badrows<-c(badvarrows,badrows)
	
	if(length(badrows)>0){
	    cat("removing bad rows\n")
	    matlist<-lapply(1:nummats, function(x){ matlist[[x]][-badrows,] })
	}
	
	numcols<-ncol(matlist[[1]])
	numgenes<-nrow(matlist[[1]])
	genenames<-row.names(matlist[[1]])
	
	#make pairwise matrix correlations
	cat("calculating global correlations\n")
	pairs<-expand.grid(1:nummats,1:nummats)
	matcors<-as.numeric(unlist(lapply(1:nrow(pairs),function(x) cor(as.vector(matlist[[pairs[x,1]]]) , as.vector(matlist[[pairs[x,2]]]) , use="complete.obs" ))))
	cormat<-matrix(nrow=nummats,ncol=nummats)
	row.names(cormat)<-legendnames
	colnames(cormat)<-legendnames
	for(i in 1:nrow(pairs)){
		cormat[pairs[i,1],pairs[i,2]]<-matcors[i]
	}
	globaltablename<-paste(matnames[1],"_globalcors.tsv",sep="")
	cat("saving global correlation table to",globaltablename,"\n")
	write.table(cormat,file=globaltablename,sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
	
	#plot pairwise global correlations as a heatmap
	globalhmname<-paste(matnames[1],"_matcors.pdf",sep="")
	cat("saving matrix cross-correlation heatmap to",globalhmname,"\n")
	pdf(file=globalhmname)
	heatmap.2(cormat,Rowv=NA,Colv=NA,dendrogram="none",trace="none",col=hcolors,breaks=brks,main="Global cross-correlations")
	
	#calculate pairwise gene correlations
	cat("calculating pairwise gene correlations\n")
	genecormat<-as.data.frame(
		lapply(1:nummats,function(y)
			unlist(
				lapply(
					1:numgenes, function(x){
						cor( matlist[[y]][x,],matlist[[1]][x,], use="complete.obs" )
					}
				)
			)
		)
		
	)
	row.names(genecormat)<-genenames
	colnames(genecormat)<-legendnames
	genecormat<-data.matrix(genecormat)
	
	genemeta<-data.frame(
	     "symbol"=unlist(lapply(strsplit(row.names(genecormat),"_"), "[", 3 ) ),
	     "id"=unlist(lapply(strsplit(row.names(genecormat),"_"), "[", 1 ) ),
	     "chrom"=unlist(lapply(strsplit(row.names(genecormat),"_"), "[", 2 ) ),
	     stringsAsFactors=FALSE
	)
	
	goodrows<-complete.cases(genecormat)
	genecormat<-genecormat[goodrows,]
	genemeta<-genemeta[goodrows,]
	
	#kmeans cluster and plot gene correlations
	cat("kmeans clustering data and plotting heatmap\n")
	k<-kmeans(genecormat,numgroups)
	groupcolors<-unlist(lapply(k$cluster,function(x) rainbow(numgroups)[x]))[order(k$cluster)]
	heatmap.2(genecormat[order(k$cluster),],trace="none",Colv=NA,Rowv=NA,dendrogram="none",RowSideColors=groupcolors,col=hcolors,breaks=brks,labRow=NA,margins=c(10,10),cexCol=1,main=paste("clustered gene-by-gene correlations with",legendnames[1]))
	cat("plotting sorted heatmap\n")
	for(i in 2:nummats){
		heatmap.2(genecormat[order(genecormat[,i]),],trace="none",Colv=NA,Rowv=NA,dendrogram="none",col=hcolors,breaks=brks,labRow=NA,main=paste("gene-by-gene correlations with ",legendnames[1],"sorted on",legendnames[sorton]))
	}
	if(hcluster==TRUE){
		cat("heirarchical clustering data and plotting heatmap\n")
		heatmap.2(genecormat,Colv="none",dendrogram="row",trace="none",col=hcolors,breaks=brks,labRow=NA,main=paste("hclustered gene-by-gene correlations with ",legendnames[1]))
	}
	
	#save file
	rbc<-rainbow(nummats-1)
	plot(density(genecormat[,2],from=-1,to=1),col=rbc[1],xlim=c(-1,1),ylim=c(0,10),xlab=paste("gene-by-gene correlation with",legendnames[1]))
	if(nummats>2){
		for(i in 3:(nummats)){
			lines(density(genecormat[,i],from=-1,to=1),col=rbc[i-1])
		}
	}
	legend("topleft",legend=legendnames[2:nummats],col=rbc,lwd=3)
	
	par(mar=c(7,5,1,1))
	boxplot(genecormat[,2:nummats],ylab=paste("gene-by-gene correlation with",legendnames[1]),las=2)
	
	dev.off()
	genecormat<-cbind(genecormat,genemeta)
	write.table(genecormat,file=paste(matnames[1],"_gene_correlations.tsv",sep=""),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
}
mat.strandtosense	<- function( topstrandmats,botstrandmats,bed){
	curbed<-read.tsv(bed)
	negrows<-which(curbed[,6]=="-")
	for(i in 1:length(topstrandmats)){
		cat("processing",topstrandmats[i],"\n")
		curplus<-read.mat(topstrandmats[i])
		curminus<-read.mat(botstrandmats[i])
		newplus<-curplus
		newminus<-curminus
		newplus[negrows,]<-curminus[negrows,]
		newminus[negrows,]<-curplus[negrows,]
		write.mat(newplus,file=gsub("\\.mat","_sense.mat",topstrandmats[i]))
		write.mat(newminus,file=gsub("\\.mat","_antisense.mat",botstrandmats[i]))
	}
}
mat.plotrows		<- function( mats,genes,geneid="symbols",sep="_",prunerows=FALSE,lspan=0,normalize=F,legendnames=basename(removeext(mats)),ylims=c(0,"auto"),plotcolors=rainbow(length(mats)),filename="orphan",start="beginning",stop="end",center="TSS",allforward=TRUE,linewidth=3,blockheight=20,grid=TRUE,X=F){
	library(tools)
	library(parallel)
	nummats<-length(mats)
	
	#read in matrices and get info
	cat("loading matrices...\n")
	matlist<-mclapply(mats,read.mat,mc.cores=detectCores())
	matcols<-unlist(lapply(matlist,ncol))
	numgenes<-length(genes)
	if(numgenes>16 & X==T){X=F;cat("more than 16 genes, forcing X=FALSE\n")}
	matnames<-basename(removeext(mats))
	genenames<-unlist(lapply(strsplit(row.names(matlist[[1]]),sep),"[",3))
	windowsizes=as.numeric( gsub("mat","",file_ext(mats) ) )
	xs<-lapply(1:length(mats),function(x) ((1:matcols[x])-(matcols[x]/2))*windowsizes[x])
	outname<-uniquefilename("mat_plotrows.pdf")
	
	if(prunerows==TRUE){
		cat("finding empty rows\n")
		badrows<-unique(unlist(lapply(1:nummats, function(x){ which(rowSums(is.na(matlist[[x]]))==matcols[x]) })))
		if(length(badrows)>0){
			cat("removing empty rows\n")
			matlist<-lapply(1:nummats,function(x){
				matlist[[x]][-badrows,]
			})
		}
	}
	
	
	#normalize data
	if(normalize==TRUE){
		cat("normalizing data\n")
		normat<-function(curmat){t(simplify2array(lapply(1:nrow(curmat),function(x) curmat[x,]/mean(curmat[x,],na.rm=TRUE))))}
		matlist<-lapply(matlist,normat)
	}
	
	#find ylims
	cm=1
	if(ylims[1]=="auto" | ylims[2]=="auto"){
		cat("gathering data to calculate ylims\n") 
		cm<-do.call(cbind,matlist)
	}
	
	if(ylims[1]=="auto"){
		cat("calculating ymins\n")
		ymins<-apply(cm,1,min,na.rm=TRUE)
		#print(ymins)
	}
	else{
		ymins<-as.numeric(rep(ylims[1],length(genenames)))
	}
	if(ylims[2]=="auto"){
		cat("calculating ymaxs\n")
		ymaxs<-1.1*apply(cm,1,max,na.rm=TRUE)
		#ok, shprint(ymaxs)
	}
	else{
		ymaxs<-as.numeric(rep(ylims[2],length(genenames)))
	}
	rm(cm)
	
	#identify genes to plot
	cat("Finding genes\n")
	if(geneid=="symbols"){
		matrows<-unlist(lapply(1:numgenes,function(x) which(genenames==genes[x])[1]))
	}
	if(geneid=="ucsc"){
		matrows<-unlist(lapply(1:numgenes,function(x) which(genenames==genes[x])[1]))
	}
	if(geneid=="rows"){
		matrows<-genes
	}
	
	if(X==F){
		pdf(file=outname)
	}
	else{
		parnum<-ceiling(sqrt(numgenes))
		par(mfrow=c(parnum,parnum))
	}
	
	#plot averages
	plot(
			0,
			type="n",
			xlim=c(min(xs[[1]]),max(xs[[1]])),
			ylim=c(mean(ymins),mean(ymaxs)),
			xlab="Distance from TSS (bp)",
			ylab="Average FPM",
			main=paste("Average score for",numgenes,"genes"),
			cex.lab=1.5,
			cex.axis=1.5,
			cex.main=1.5,
			cex.sub=1.5
	)
		
	for(m in 1:length(mats)){
			x<-1:matcols[m]
			y<-colMeans(matlist[[m]][matrows,],na.rm=TRUE)
			lines(
				xs[[m]],
				y,
				col=plotcolors[m],
				lwd=3
			)
	}
		
	legend("topright",legend=legendnames,col=plotcolors,lwd=3)
	
	
	
	for(i in 1:numgenes){
		plot(
			0,
			type="n",
			xlim=c(min(xs[[1]]),max(xs[[1]])),
			ylim=c(ymins[matrows[i]],ymaxs[matrows[i]]),
			xlab="Distance from TSS (bp)",
			ylab="Normalized FPM",
			main=genenames[matrows[i]],
			cex.lab=1.5,
			cex.axis=1.5,
			cex.main=1.5,
			cex.sub=1.5
		)
		
		for(m in 1:length(mats)){
			if(lspan!=0){
				x<-1:matcols[m]
				y<-matlist[[m]][matrows[i],]
				y[-which(is.na(y))]<-loess(y~x,span=lspan)$fitted
			}
			else{
				y<-matlist[[m]][matrows[i],]
			}
			lines(
				xs[[m]],
				y,
				col=plotcolors[m],
				lwd=3
			)
		}
		
		legend("topright",legend=legendnames,col=plotcolors,lwd=3)
	}
	cat("pdf saved to",outname,"\n")
	
	if(X==F){dev.off()}
}
mat.colcors		<- function( matfile1,matfile2,center="TSS",continue=FALSE,color="black"){
	
	mat1<-data.matrix(read.table(paste(matfile1),stringsAsFactors=FALSE,row.names=NULL))
	mat2<-data.matrix(read.table(paste(matfile2),stringsAsFactors=FALSE,row.names=NULL))
	correlations<-as.numeric(vector(length=ncol(mat1)))
	for(i in 1:ncol(mat1)){
		correlations[i]<-cor(mat1[,i],mat2[,i],use="p")
	}

	if(continue==FALSE){plot(((1:ncol(mat1))-(ncol(mat1)/2)),correlations,xlab=paste("Distance from",center),ylab="Correlation",type="l",col=color,main=paste("Correlation between",removeext(matfile1),"and",removeext(matfile2),"along",center),ylim=c(0,1))}
	else{lines(((1:ncol(mat1))-(ncol(mat1)/2)),correlations,xlab=paste("Distance from",center),ylab="Correlation",col=color,main=paste("Correlation between",removeext(matfile1),"and",removeext(matfile2),"along",center),ylim=c(0,1))}
}
mat.plotaverages		<- function( matrixlist,centername="TSS",plottype="l",ylims=c("auto","auto"),X=TRUE,prefix="unnamed",plotcolors=rainbow(length(matrixlist)),legendnames=basename(removeext(matrixlist))){
	library(tools)
	library(parallel)
	cores=detectCores()-1
	nummats<-length(matrixlist)
	matlist<-mclapply(matrixlist,read.mat,mc.cores=cores)
	mataverages<-mclapply(matlist,colMeans,na.rm=TRUE,mc.cores=cores)
	ymax<-max(unlist(mataverages),na.rm=TRUE)
	ymin<-min(unlist(mataverages),na.rm=TRUE)
	if(ylims[1]=="auto"){ylims[1]=ymin}
	if(ylims[2]=="auto"){ylims[2]=ymax}
	ylims<-as.numeric(ylims)
	windowsizes<-as.numeric(gsub("mat","",file_ext(matrixlist)))
	matcols<-unlist(lapply(matlist,ncol))
	xs<-lapply(1:nummats,function(x) ((1:matcols[x])-(matcols[x]/2))*windowsizes[x])
	if(X==FALSE){pdf(file=paste(prefix,"_averages.pdf",sep=""))}
	plot(
        	0,
		type="n",
    		xlim=c(min(xs[[1]]),max(xs[[1]])),
		ylim=ylims,
    		xlab=paste("Distance from",centername,"(bp)"),
    		ylab="Average score",
    		cex=0.5
		#cex.lab=1.5,
    		#cex.axis=1.5,
    		#cex.main=1.5,
    		#cex.sub=1.5
    	)
	for(k in 1:nummats){
		lines(
		      xs[[k]],
		      mataverages[[k]],
		      col=plotcolors[k],
		      lwd=1,
		      type=plottype
		)
    	}
	legend("topleft",legend=legendnames, col=plotcolors, lwd=3)
	if(X==FALSE){dev.off()}
}
mat.plotaverages2	<- function( matfilelist,suffixes,legendnames,plotcolors=rainbow(length(matfilelist)), windowsize=10, matcols=100, xname="Distance from Feature", yname="Average Score"){
	numlists<-length(matfilelist)
	if(length(suffixes==1)){suffixes<-rep(suffixes,numlists)}
	prefixes<-lapply(1:numlists, function(x){
		remove.suffix(basename(matfilelist[[x]]),suffixes[x])
	})
	feats<-unique(unlist(prefixes))
	matches<-lapply(1:numlists,function(x) match(feats,prefixes[[x]]) )
	pdfname<-paste("mataverages_",length(feats),"centers_",numlists,"scores.pdf",sep="")
	pdf(file=pdfname)
	xaxis<-windowsize*(1:matcols) - (windowsize*matcols)/2
	for(i in 1:length(feats)){
		cat("plotting feature #",i,"/",length(feats),":",feats[i],"\n")
		refs<-unlist(lapply(matches,"[",i))
		nomats<-which(is.na(refs))
		goodmats<-which(1:numlists %ni% nomats)
		curmats<-lapply(goodmats,function(x) read.mat(matfilelist[[x]][matches[[x]][i]]))
		counts<-unlist(lapply(curmats,nrow))
		colmeans<-lapply(curmats,colMeans,na.rm=TRUE)
		all<-unlist(colmeans)
		ylims<-c(min(all),max(all))
		plot(0,type="n",xlim=c(min(xaxis),max(xaxis)),ylim=ylims,main=feats[i],xlab=xname,ylab=yname)
		n<-rep(0,numlists)
		n[goodmats]<-counts
		for(j in 1:numlists){
			if(j %in% goodmats == TRUE){ 
				lines(xaxis,colmeans[[which(goodmats==j)]],col=plotcolors[j],lwd=3)
			}
			
		}
		legend("topleft",legend=paste(legendnames,n),col=plotcolors,lwd=3)
	}
	dev.off()
	cat("pdf saved to",pdfname,"\n")
}
mat.plotsamplerows	<- function( matname,samples=5,ylims=c(0,30)){
	mat<-read.mat(matname)
	rows<-sample(1:nrow(mat),samples)
	plotcolors=rainbow(samples)
	plot(1:ncol(mat),mat[rows[1],],col=plotcolors[1],type="l",ylim=ylims)
	for(i in 2:samples){
		lines(1:ncol(mat),mat[rows[i],],col=plotcolors[i])
	}
	legend("topright",legend=row.names(mat)[rows],col=plotcolors,lwd=2)
}
mat.loess		<- function( mats, lspan=0.05 ){
	library(parallel)
	matnames<-basename(removeext(mats))
	print(matnames)
	for(i in 1:length(mats)){
		cat("reading mat\n")
		m<-read.mat(mats[i])
		y<-1:ncol(m)
		ml<-data.matrix(t(as.data.frame(mclapply(1:nrow(m),function(x){  m[x,-which(is.na(m[x,]))]<-loess(m[x,]~y,span=lspan)$fitted ; m[x,]   },mc.cores=detectCores()))))
		cat("smoothing\n")
		row.names(ml)<-row.names(m)
		cat("saving mat\n")
		write.mat(ml,file=paste(matnames[i],"_loess",gsub("\\.","-",as.character(lspan)),".mat10",sep=""))
	}
}
mat.transform		<- function( mats, from="log2", to="ratio"){}
mat.qnorm		<- function( mats, normalizeto=1 ){
	library(tools)
	library(parallel)
	matnames<-basename(removeext(mats))
	exts<-file_ext(mats)
	cat("loading matrices\n")
	matlist<-mclapply(mats,read.mat,mc.cores=detectCores())
	nummats<-length(mats)
	cat("normalizing matrices\n")
	matlist<-mclapply(1:nummats,function(x) {matlist[[x]][order(matlist[[x]])]<-matlist[[normalizeto]][order(matlist[[normalizeto]])];matlist[[x]]} )
	cat("saving matrices\n")
	outnames<-paste(matnames,"_qnorm.",exts, sep="") 
	mclapply(1:nummats, function(x) write.mat(matlist[[x]],file=outnames[x]) ,mc.cores=detectCores() )
	return(outnames)
}
mat.hist			<- function( mats, xlims=c("auto","auto"), plotcolors=rainbow(length(mats)), legendnames=basename(removeext(mats)) ){
	library(parallel)
	nummats<-length(mats)
	
	cat("reading in matrices\n")
	matlist<-mclapply( mats, read.mat, mc.cores=detectCores() )
	
	cm=1
	if(xlims[1]=="auto" | xlims[2]=="auto"){
		cat("calculating xlims\n")
		cm<-do.call(cbind,matlist)
		q<-quantile(cm,probs=c(0.05,0.95),na.rm=T)
	}
	if(xlims[1]=="auto"){xlims[1]=q[1]}
	if(xlims[2]=="auto"){xlims[2]=q[2]}
	xlims<-as.numeric(xlims)
	print(xlims)
	
	cat("calculating score distributions\n")
	densitylist<-lapply( 1:nummats, function(x) density( matlist[[x]], from=xlims[1], to=xlims[2], na.rm=T))
	ymax<-max(unlist(lapply( 1:nummats, function(x) density( matlist[[x]], from=xlims[1], to=xlims[2],na.rm=T)$y  )),na.rm=TRUE)
	print(ymax)
	cat("plotting score distributions\n")
	plot(0,type="n",ylim=c(0,ymax), xlim=xlims, xlab="Fragment Size", ylab="Density")
	for(i in 1:nummats){
		lines(densitylist[[i]],col=plotcolors[i],lwd=2)
	}
	legend("topright",legend=legendnames, col=plotcolors, lwd=3)
}
mat.window		<- function( matfile, operation="mean", windowsize=5){
    library(tools)
    if(ceiling(windowsize/2)==floor(windowsize/2)){stop("windowsize must be an odd number")}
    cat("reading in matrix\n")
    mat<-read.mat(matfile)
    
    cat("gathering info\n")
    flank<-(windowsize-1)/2
    numcols<-ncol(mat)
    matrownames<-row.names(mat)
    
    print(numcols)
    
    cat("windowing matrix\n")
    mat<-simplify2array(lapply(1:numcols,function(x){
	if(x-flank<1){l=1}
	else{l=x-flank}
	if(x+flank>numcols){r=numcols}
	else{r=x+flank}
	rowMeans(mat[,l:r],na.rm=T)
    }))
    
    cat("wrapping up\n")
    row.names(mat)<-matrownames
    ext<-file_ext(matfile)
    outname<-paste(basename(removeext(matfile)),"_",operation,"win",windowsize,".",ext,sep="")
    write.mat(mat,file=outname)
    return(outname)
}
mat.getscores		<- function( matfile, refmatfile ){
	mat1<-read.mat(matfile)
	mat2<-read.mat(refmatfile)
	mat1genes<-unlist(lapply(strsplit(row.names(mat1),"_"),"[",1))
	mat1genes<-unlist(lapply(strsplit(row.names(mat1),"-"),"[",1))
	mat2genes<-unlist(lapply(strsplit(row.names(mat2),"_"),"[",1))
	refrows<-match(mat1genes,mat2genes)
	badrows<-which(is.na(refrows))
}
mat.heatmap3		<- function( mats=NULL , matnames=NULL , sorton=1 , sort.methods="none" , gene.list = NULL , sort.ranges = NA , sort.breaks = NA , numgroups=3 , crossmat=NULL , matcolors=NULL , groupcolors=NULL , heatmap.lims="5%,95%" , heatmap.colors="black,white" , fragmats=NULL , fragmatnames=NULL , fragrange=c(0,200) , vplot.colors="black,blue,yellow,red" , vplot.lims="0,95%" , forcescore=TRUE , cores="max" ){
	
	
	# TO DO
	# ##############
	# check sorting for expected sorting methods or existing file
	# check that mats is a character vector and all files exist
	# check that all mats have identical row names and nrow
	# check that all colors are R colors

	library(gtools)
	library(tools)
	library(parallel)
	pwd<-getwd()
	if(cores=="max"){cores=detectCores()-1}
	if(is.null(mats) & is.null(fragmats)){stop("nothing to plot")}
	

	valid.sort.methods=c("none","kmeans","chrom","max","min","mean","median","sd","maxloc","minloc","pregrouped")
	if(sum(sort.methods %in% valid.sort.methods) != length(sort.methods)){stop("unrecognized sort.methods. valid methods include none,kmeans,chrom,max,min,mean,median,sd,maxloc,minloc")}

	#CHECK SCOREMAT PARAMETERS AND DEFINE SCOREMAT SETTINGS
	if(is.null(mats) == FALSE){
		
		#gather info and check paramters
		nummats<-length(mats)
		winsizes=as.numeric(gsub("mat","",file_ext(mats)))
		if(sum(is.na(winsizes))>0){stop("cannot find window sizes for one or more matrices. window size (in bp) should be in the matrix file extension after 'mat' (for example, matrixname.mat10 for a window size of 10 bp")}
		# #### check mat names for uniqueness and length

		# set the names of matrices using file names if they aren't explicitly defined
		if(is.null(matnames)){ matnames<-basename(removeext(mats)) }
		
		#set unspecified settings to defaults
		if(length(heatmap.lims)==1){heatmap.lims<-rep(heatmap.lims,nummats)}
		

		if(length(heatmap.colors)==1){heatmap.colors<-rep(heatmap.colors,nummats)}
		scolramps<-lapply(heatmap.colors, function(x) colorRampPalette(unlist(strsplit(x,","))))
		
	}	
	
	#CHECK FRAGMAT PARAMETERS AND DEFINE FRAGMAT SETTINGS
	if(is.null(fragmats) == FALSE){
		numfragmats<-length(fragmats)
		fragwinsizes=as.numeric(gsub("fmat","",file_ext(fragmats)))
		if(sum(is.na(fragwinsizes))>0){stop("cannot find window sizes for one or more matrices")}
		if(is.null(fragmatnames)){ fragmatnames<-basename(removeext(fragmats)) }
		if(length(vplot.lims)==1){vplot.lims<-rep(vplot.lims,numfragmats)}
		

		if(length(vplot.colors)==1){vplot.colors<-rep(vplot.colors,numfragmats)}
		vcolramps<-lapply(vplot.colors, function(x) colorRampPalette(unlist(strsplit(x,","))) )
		
		vplotnames1<-vector(mode="character")
		vplotnames2<-vplotnames1
		

	}
	
	
	
	
	#load sort table if defined, otherwise sort based on mats defined in sorton
	if(mode(sort.methods)=="file"){
		# check if file exists and if sort.methods is length 1
		cat("loading sort table\n")
		sorttable<-read.table(sorting,quote=F,row.names=1,header=T,sep="\t",stringsAsFactors=F)
		sortsettings<-data.frame("V1"=colnames(sorttable)[1:(ncol(sorttable)-2)],"V2"="presort","V3"="presort")
		numgroups<-unlist(lapply(sorttable,length)) #THIS PROBABLY WRONG
		numsorts<-nrow(sortsettings)
		
		#make directory for saving files
		if(file.exists("heatmaps")==FALSE){system(command="mkdir heatmaps")}
		dname<-paste("heatmaps/",basename(removeext(r)),sep="")
		#dname<-uniquefilename(dname)
		system(paste("mkdir",dname))
		write.table(sorttable,file=paste(dname,"/",basename(dname),".sort",sep=""),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
	} else{
		
		# read matrix for sorting and get dimensions
		sortmat<-read.mat(mats[sorton])
		matcols<-ncol(sortmat)
		matrows<-nrow(sortmat)
		
		# tabulate sort criteria
		sortsettings<-data.frame(
			"METHOD"=sort.methods,
			"START"=as.numeric(unlist(lapply(strsplit(sort.ranges,","),"[",1))),
			"END"=as.numeric(unlist(lapply(strsplit(sort.ranges,","),"[",2))),
			stringsAsFactors=FALSE
		)

		# if only 1 sort method, repeat method so sorting script runs OK
		if(nrow(sortsettings)==1 & sortsettings[1,1] != "kmeans" & sortsettings[1,1] != "chrom"){
			sortsettings[2,1]<-sortsettings[1,1]
		}

		# if sort range not defined, specify to use entire matrix
		sortsettings[is.na(sortsettings[,2]),2]=as.numeric((-1*winsizes[sorton])*(matcols-matcols/2))
		sortsettings[is.na(sortsettings[,3]),3]=as.numeric(winsizes[sorton]*(matcols-matcols/2))
		
		# determine how many sort criteria are defined, and match length of vector defining how many groups to create for each sort criterion
		numsorts<-nrow(sortsettings)
		length(numgroups)<-numsorts
		numgroups[is.na(numgroups)]<-1
		
		#print sort settings
		cat("\n####################\n")
		cat("SORT CRITERIA\n")
		print(sortsettings)
		cat("####################\n\n")
		
		#convert bp to column #s in sort settings
		sortsettings[,2]<-sortsettings[,2]/winsizes[sorton] + matcols/2 +1
		sortsettings[,3]<-sortsettings[,3]/winsizes[sorton] + matcols/2
		
		#if(normalize==TRUE){
		#	cat("normalizing sort matrix\n")
		#	sortmat<-t(simplify2array(mclapply(1:matrows[r],function(x) sortmat[x,]/mean(sortmat[x,],na.rm=TRUE),mc.cores=cores)))*mean(sortmat,na.rm=TRUE)
		#	row.names(sortmat)<-matlist[[r]]
		#}


		#sorttable<-data.frame("V1"=rep(1,matrows),stringsAsFactors=F)
		sorttable<-data.frame(matrix(NA,nrow=matrows,ncol=numsorts),stringsAsFactors=F)
		grouptable<-sorttable

		for(s in 1:numsorts){

			#sortcol<-vector(mode="numeric",length=matrows)

			if(s==1){numprevgroups=1} else{ numprevgroups<-length(unique(grouptable[,s-1])) }
			
			for(p in 1:numprevgroups){
				
				if(s==1){grouprows=1:matrows} else{grouprows<-which(grouptable[,s-1]==p)}
				
				#make and trim sort matrix
				tmpsortmat<-sortmat
				tmpsortmat<-tmpsortmat[grouprows,sortsettings[s,2]:sortsettings[s,3]]
				submatrows<-nrow(tmpsortmat)

				#convert infinite values to NA
				tmpsortmat[is.infinite(tmpsortmat)]<-NA
				
				#perform calculation for sorting
				if(sortsettings[s,1] == "sd"){
					cat("calculating standard deviation\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,sd, na.rm=TRUE)))
				}
				if(sortsettings[s,1] == "median"){
					cat("calculating medians\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,median, na.rm=TRUE)))
				}
				if(sortsettings[s,1] == "max"){
					cat("calculating maxs\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,max, na.rm=TRUE)))
				}
				if(sortsettings[s,1] == "min"){
					cat("calculating mins\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,min, na.rm=TRUE)))
				}
				if(sortsettings[s,1] == "mean"){
					cat("calculating means\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,mean, na.rm=TRUE)))
				}
				if(sortsettings[s,1] == "maxloc"){
					cat("calculating max locations\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,which.max)))
				}
				if(sortsettings[s,1] == "minloc"){
					cat("calculating min locations\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,which.min)))
				}
				if(sortsettings[s,1] == "chrom"){
					cat("sorting by chromosome\n")
					chromlist<-unlist(lapply(strsplit(row.names(sortmat),";"),"[",2))
					chroms<-unique(chromlist)
					chroms<-mixedsort(chroms)
					chromnums<-1:length(chroms)
					numgroups[s]<-length(chroms)
					sortvals<-match(chromlist,chroms)
				}
				if(sortsettings[s,1] == "kmeans"){
					cat("kmeans-clustering data (k=",numgroups[s],")\n",sep="")
					tmpsortmat2<-tmpsortmat
					tmpsortmat2[is.na(tmpsortmat2)]<-0
					sortvals<-kmeans(tmpsortmat2,numgroups[s])$cluster
					rm(tmpsortmat2)
				}
				if(sortsettings[s,1] == "none"){
					cat("not sorting\n",sep="")
					sortvals<-1:nrow(tmpsortmat)
				}
				if(sortsettings[s,1] == "pregrouped"){
					cat("not sorting\n",sep="")
					sortvals<-tmpsortmat[,1]
					numgroups[s]<-length(unique(sortvals))
				}
				

				#assign sortvals to sorttable
				sorttable[grouprows,s] <- sortvals
				
				# assigned pregrouped data to group table
				if(sortsettings[s,1] == "kmeans" | sortsettings[s,1] == "chrom" | sortsettings[s,1] == "pregrouped"){
					grouptable[grouprows,s] <- sortvals
				}

				# if not already pregrouped, group sort values into specified # of groups
				if(sortsettings[s,1] != "kmeans" & sortsettings[s,1] != "chrom" | sortsettings[s,1] != "pregrouped"){
					
					# #### what if numgroups[s]==1 ???????

					# if breaks not defined for a given sorting criteria
					if(is.na(sort.breaks[s])){
						
						
						# find number of genes in an equally-split group
						rowspergroup<-ceiling(submatrows/numgroups[s])

						# define breaks based on number of genes per group
						brks<-c(seq(1,submatrows,rowspergroup),submatrows)

						#assign rows into equally-sized groups
						#grouptable[grouprows[order(sortvals)],s]<-cut((1:submatrows)[order(sortvals)],brks,labels=F,include.lowest=T)
						grouptable[grouprows,s]<-cut((1:submatrows)[order(order(sortvals))],brks,labels=F,include.lowest=T)

					} else{
						
						cat("grouping genes into",numgroups[s],"groups based on breaks\n")
					
						# #set numgroups to length(brks)+1
						numgroups[s]<-length(brks)+1
						
						brks<-unlist(strsplit(sort.breaks[s],","))
						
						# if percentiles defined for breaks, calculate break values
						if(grepl("%",brks)){
							
							# calculate quantile from %s in sort.breaks
							brks[grep("%",brks)]<-quantile(sortvals,probs=as.numeric(gsub("%","",brks[grep("%",brks)])))
						}
						brks<-c(-Inf,brks,Inf)
						
						# group genes by defined breaks
						grouptable[grouprows,s] <- cut(sortvals,brks,include.lowest=T,labels=F)
					}
				}
			}
		}

		# create table of colors
		if(is.null(groupcolors)){ groupcolors<-rainbow(max(numgroups)) }
		getcolor <- function(c){groupcolors[c]}
		color.table <- as.data.frame ( lapply ( grouptable , getcolor ) )

		# sort table
		heatmap.order <- order(do.call(order,as.data.frame(sorttable)))

		# name columns of tables
		colnames(grouptable)<-paste("GROUPING",1:numsorts,"-",sortsettings[,1],sortsettings[,3],sortsettings[,3],sep=";")
		colnames(sorttable)<-paste("SORTVALS",1:numsorts,"-",sortsettings[,1],sortsettings[,3],sortsettings[,3],sep=";")
		colnames(color.table)<-paste("COLORS",1:numsorts,"-",sortsettings[,1],sortsettings[,3],sortsettings[,3],sep=";")
		
		#make a table to store order, group, and other info
		heatmap.table <- data.frame(cbind(
			row.names=row.names(sortmat),
			heatmap.order,
			sorttable,
			grouptable,
			color.table
		))
		
		# ###### !! NEED FIXING AFTER MAT.MAKE IS CHANGED ####################
		#sorttable$symbol<-unlist(lapply(strsplit(row.names(sorttable),";"), "[", 3 ) )
		#sorttable$id<-unlist(lapply(strsplit(row.names(sorttable),";"), "[", 1 ) )
		#sorttable$chromosome<-unlist(lapply(strsplit(row.names(sorttable),";"), "[", 2 ) )		
		
		
		#make directory for saving files
		dir.create("heatmaps", showWarnings = FALSE)
		dname<-paste("heatmaps/",basename(matnames[sorton]),"_",numgroups[1],"g_",sortsettings[1,1],"_",sortsettings[1,2],"-",sortsettings[1,3],sep="",collapse="")
		#dname<-uniquefilename(dname)
		system(paste("mkdir",dname))
		write.table(heatmap.table,file=paste(dname,"/",basename(matnames[sorton]),".sort",sep=""),sep="\t",row.names=TRUE,quote=FALSE,col.names=TRUE)
	}
	
	
	if(is.null(mats)==FALSE){

		# ##### get order of preloaded sort table
		
		# get number of rows of sort table if preloaded
		matrows<-nrow(sorttable)

		# order the table of group info
		heatmap.order<-heatmap.table$heatmap.order
		table.starts<-1:3*(ncol(heatmap.table)-1)/3
		grouptable<-as.matrix(heatmap.table[,table.starts[2]:(table.starts[2]+1)])
		grouptable<-grouptable[order(heatmap.order),]
		color.table<-as.matrix(heatmap.table[,table.starts[3]:(table.starts[3]+1)])
		color.table<-color.table[order(heatmap.order),]
		#grouptable<-as.matrix(heatmap.table[,table.starts[2]:(table.starts[2]+1)])[order(heatmap.order,na.last=F),]
		#color.table<-as.matrix(heatmap.table[,table.starts[3]:(table.starts[3]+1)])[order(heatmap.order,na.last=F),]
		#heatmap.table<-heatmap.table[order(heatmap.order,na.last=F),]
		
		# identify number of rows in each group for each sorting
		clusternums<-lapply(1:numsorts,function(h){ unlist(lapply(1:numgroups[h],function(x) length(which(grouptable[,h]==x)))) } )
		grouprownumbers<-lapply(1:numgroups[1], function(x) which(grouptable[,1]==x) )
		
		# set legend names and colors
		if(numgroups[1]>1){
			legendnames<-paste(c(1:numgroups[1],"all"),", n= ",c(clusternums,matrows),sep="")
			
			#FIX THIS ######################################
			if(sortsettings[1,1]=="chrom"){
				legendnames<-paste(c(chroms,"all"),", n= ",c(clusternums,matrows),sep="")
			}
			
			legendcolors<-c(groupcolors,"gray50")

		} else{
			legendnames<-paste("all, n =",matrows)
			legendcolors<-"black"
		}
	}

	#process each frag matrix
	if(is.null(fragmats) == FALSE){
		for(v in 1:numfragmats){
			

			fragmat<-read.fmat(fragmats)
			
			fmatcols<-ncol(fragmat)

			x<-((1:fmatcols)-(fmatcols/2))*fragwinsizes[v]

			fragsizes<-lapply(1:fmatcols,function(x) as.numeric(na.omit(as.numeric(unlist(strsplit(paste(fragmat[,x],collapse=","),","))))))
			
			fragsizes<-lapply(1:fmatcols,function(x) fragsizes[[x]][which(fragsizes[[x]] >= fragrange[1] & fragsizes[[x]] <= fragrange[2])])

			vmat<-matrix(0,ncol=fmatcols,nrow=fragrange[2]-fragrange[1])
			
			for(h in 1:fmatcols){
				vmat[,h]<-hist(fragsizes[[h]],breaks=fragrange[2]:fragrange[1],plot=F)$counts
			}
				
			#set vmat score limits
			vbot<-strsplit(vplot.lims[v],",")[[1]][1]
			vtop<-strsplit(vplot.lims[v],",")[[1]][2]
			if(grepl("%",vbot)){ vbot<-quantile(vmat , probs=as.numeric(gsub("%","",vbot)) /100,na.rm=T) }
			if(grepl("%",vtop)){ vtop<-quantile(vmat , probs=as.numeric(gsub("%","",vtop)) /100,na.rm=T) }
			vbot<-as.numeric(vbot)
			vtop<-as.numeric(vtop)
			if(vtop==vbot){vtop=vtop+1}
			
			# if(rpm==TRUE){
			# 	vscalar<-1000000/as.numeric(gsub("#","",readLines(pipe(paste("head -n 1",fragmats[v])))))
			# 	if(is.na(vscalar)==FALSE){vmat<-vmat*vscalar
			# 	} else{cat("WARNING: FRAGMENT COUNT NOT FOUND IN FRAGMAT, NOT NORMALIZING\n")}
			# }

			#vmat<-vmat[nrow(vmat):1,]
			
			cat("drawing vplot\n")
			
			brks=seq(vbot,vtop,by=(vtop-vbot)/100)

			

			# draw heatmap of all
			vplotname<-paste(pwd,"/",dname,"/","vplot_all_",fragmatnames[v],".png",sep="")
			png(width=1000,height=1000,file=vplotname)
			

			
			layout(matrix(c(4,2,1,3),nrow=2),heights=c(1,4),widths=c(4,1))
			#layout(matrix(1:2,nrow=2),heights=c(1,4))
			par(oma=c(2,2,2,2))
			par(mar=c(3,3,0,0))
			image(matrix(brks),col=vcolramps[[v]](100),breaks=brks,axes=F)
			axis(side=1,at=c(0,1),labels=c(brks[1],brks[length(brks)]))
			mtext("bin counts",side=1,line=2)
			
			par(mar=c(3,3,0,0))
			image(t(vmat),breaks=brks,col=vcolramps[[v]](100),xaxt='n',yaxt='n')
			axis(side=1,at=c(0,0.5,1),labels=c(min(x),0,max(x)))
			axis(side=2,at=(0:4)/4,labels=fragrange[1]+0:4*((fragrange[2]-fragrange[1])/4))
			mtext("Distance from feature",side=1,line=2)
			mtext("Fragment size (bp)",side=2,line=2)
			
			dx<-apply(vmat,2,mean,na.rm=T)
			dy<-apply(vmat,1,mean,na.rm=T)
			par(mar=c(3,0,0,0))
			plot(dy,1:(fragrange[2]-fragrange[1]),type="l",yaxs='i',axes=F,lwd=4)
			par(mar=c(0,3,1,0))
			plot(1:fmatcols,dx,type="l",xaxs='i',axes=F,lwd=4,main=fragmatnames[v])
			

			dev.off()
			
			if(is.null(mats)==FALSE){

				fragmat<-fragmat[order(heatmap.order),]

				for(q in 1:numgroups[1]){
					

					vplotname<-paste(pwd,"/",dname,"/","vplot_group",q,"_",fmatname,".png",sep="")
					vplotname<-paste(pwd,"/",dname,"/","vplot_group",q,"_",fragmatnames[v],"_pergene.png",sep="")
					

					groupmat<-fragmat[grouprownumbers[[q]],]
					
					fragsizes<-lapply(1:fmatcols,function(x) as.numeric(na.omit(as.numeric(unlist(strsplit(paste(groupmat[,x],collapse=","),","))))))
			
					fragsizes<-lapply(1:fmatcols,function(x) fragsizes[[x]][which(fragsizes[[x]] >= fragrange[1] & fragsizes[[x]] <= fragrange[2])])

					vmat<-matrix(0,ncol=fmatcols[v],nrow=fragrange[2]-fragrange[1])
					
					for(h in 1:fmatcols[v]){
						vmat[,h]<-hist(fragsizes[[h]],breaks=fragrange[2]:fragrange[1],plot=F)$counts
					}

					# if(rpm==TRUE){
					# 	vscalar<-1000000/as.numeric(gsub("#","",readLines(pipe(paste("head -n 1",fragmats[v])))))
					# 	if(is.na(vscalar)==FALSE){vmat<-vmat*vscalar
					# 	} else{cat("WARNING: FRAGMENT COUNT NOT FOUND IN FRAGMAT, NOT NORMALIZING\n")}
					# }

					vmat<-vmat[nrow(vmat):1,]

					pergenescale<-nrow(fragmat)/nrow(groupmat)
					
					brks=seq(vbot,vtop,by=(vtop-vbot)/100)
					
					print(brks)
					print(length(brks))

					png(width=1000,height=1000,file=vplotname)
					heatmap.2(vmat,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=brks,labRow=NA,labCol=NA,col=vcolramps[[v]](100),main=fragmatnames[v],cex.main=5)
					dev.off()
					png(width=1000,height=1000,file=vplotname)
					heatmap.2(vmat*pergenescale,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=brks,labRow=NA,labCol=NA,col=vcolramps[[v]](100),main=fragmatnames[v],cex.main=5)
					dev.off()
				}

			}

		}
		montagecols<-numgroups[1]+1
		print(vplotnames1)
		cat("\n\n\n")
		print(vplotnames2)
		#cat(paste("montage -tile +",montagecols,"+",numfmats," ",paste(vplotnames1,collapse=" ")," ",pwd,"/",dname,"/vplot_montage_abs.png",sep=""))
		cat("\n\n\n\n")
		#cat(paste("montage -tile +",montagecols,"+",numfmats," ",paste(vplotnames2,collapse=" ")," ",pwd,"/",dname,"/vplot_montage_per.png",sep=""))
		cat("\n\n\n\n")
		
		
		system(paste("montage -geometry +2+2 -tile ",montagecols,"x",numfmats," ",paste(vplotnames1,collapse=" ", sep=" ")," ",pwd,"/",dname,"/vplot_montage_abs.png",sep=""))
		system(paste("montage -geometry +2+2 -tile ",montagecols,"x",numfmats," ",paste(vplotnames2,collapse=" ", sep=" ")," ",pwd,"/",dname,"/vplot_montage_per.png",sep=""))
	}
	
	
	
	
	# ##############################################################################
	# TEMPORARY IMPLEMENTATION OF CROSS-MATRIX GROUP AVERAGES
	#grouplist=lapply(1:nummats, function(x){lapply(1:numgroups[1],function(x){})})
	#q=1
	# ##############################################################################
	
	
	
	
	#process each matrix
	hmnames<-paste(pwd,"/",dname,"/","heatmap_",basename(matnames),".png",sep="")
	
	#for(l in 1:nummats){
	mclapply(1:nummats,function(l){
		
		# load matrix and get dimensions
		cat("reading in matrices\n")
		curmat<-read.mat(mats[l])
		matcols<-ncol(curmat)
		matrows<-nrow(curmat)
		
		cat("1\n")
		
		#sort matrix
		curmat<-curmat[order(heatmap.order),]

		# define heatmap color gradient function
		cat("2\n")
		
		# define x axis based on winsize and ncol assuming non-meta
		x<-((1:matcols)-(matcols/2))*winsizes[l]
		

		# set score limits for heatmap colors
		heatmap.mins<-unlist(lapply(strsplit(heatmap.lims,","),"[",1))
		heatmap.maxs<-unlist(lapply(strsplit(heatmap.lims,","),"[",2))
		if(grepl("%",heatmap.mins[l])){ heatmap.mins[l]<-quantile(curmat,probs=as.numeric(gsub("%","",heatmap.mins[l])) /100,na.rm=T) }
		if(grepl("%",heatmap.maxs[l])){ heatmap.maxs[l]<-quantile(curmat,probs=as.numeric(gsub("%","",heatmap.maxs[l])) /100,na.rm=T) }
		curmin<-as.numeric(heatmap.mins[l])
		curmax<-as.numeric(heatmap.maxs[l])
		
		if(curmin==curmax){curmax<-max(curmat,na.rm=TRUE)}
		if(curmin==curmax){curmax<-curmin+1}

		cat("3\n")
		

		curmat[is.infinite(curmat)]<-NA



		#calculate aggregate profiles
		groupmeans<-lapply(1:numgroups[1], function(y){
				colMeans(as.matrix(curmat[grouprownumbers[[y]],]),na.rm=TRUE)
		})
		groupmeans[[numgroups[1]+1]]<-colMeans(curmat,na.rm=TRUE)
		groupmeanmaxs<-unlist(lapply(groupmeans,max,na.rm=TRUE))
		groupmeanmins<-unlist(lapply(groupmeans,min,na.rm=TRUE))
		
		cat("4\n")
		
		# ##############################################################################
		# TEMPORARY IMPLEMENTATION OF CROSS-MATRIX GROUP AVERAGES
		#if(is.null(crossmat)==FALSE & l %in% crossmat == TRUE){
		#	grouplist[[q]]=groupmeans
		#	q=q+1
		#}
		# ##############################################################################
		
		#make NAs 0 for aesthetics
		if(forcescore==TRUE){
			curmat[is.na(curmat)]<-0
			curmat[is.nan(curmat)]<-0
			curmat[is.infinite(curmat)]<-0
		}
		
		
		
		cat("drawing heatmap (",basename(removeext(mats[l])),")\n")
		png(file=hmnames[l],height=5000,width=1000)
		
		# define color breaks
		brks<-c(seq(curmin,curmax,by=(curmax-curmin)/100),curmax)
		
		
		# copy matrix into new object
		mat<-curmat

		# create layout matrix
		m<-matrix(1:3,nrow=3)

		# set out-of-boundary scores to boundaries
		mat[which(mat>curmax)]<-curmax
		mat[which(mat<curmin)]<-curmin

		# create image layout and set margins
		layout(m,heights=c(1,20,5),widths=1)
		par(oma=c(10,0,10,0))
		par(mar=c(10,15,5,15))

		# draw color scale
		image(matrix(brks),col=scolramps[[l]](101),breaks=brks,axes=F)
		axis(side=1,at=c(0,1),labels=c(curmin,curmax),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
		
		# title image
		mtext(matnames[l],side=3,cex=8,outer=T)
		
		# set margins and draw heatmap
		par(mar=c(5,15,5,15),xpd=TRUE)
		image(t(mat),breaks=brks,col=scolramps[[l]](101),axes=FALSE)
		
		# create group color bars next to heatmap
		for(g in 1:numsorts){
			for(i in 1:numgroups[g]){
				linevec<-grouptable[,g]
				linevec[which(grouptable[,g] != i)]<-NA
				linecoords=(1:matrows)/matrows
				linecoords[which(is.na(linevec))]<-NA
				lines(rep(-0.01-g*0.05,matrows),linecoords[matrows:1],lwd=40,col=groupcolors[i],lend=1)
			}
		}
		axis(side=1,at=c(0,0.5,1),labels=F,lwd=10,padj=1,line=1,tcl=-3)
		
		# draw average plot
		par(mar=c(20,15,0,15))
		plot(0,type="n",xlim=c(min(x),max(x)),ylim=c(min(groupmeanmins),max(groupmeanmaxs)),xaxs="i",axes=FALSE)
		for(i in 1:(numgroups[1]+1)){
			lines(x,groupmeans[[i]],lwd=10,col=legendcolors[i])
		}
		axis(side=1,at=c(min(x),0,max(x)),labels=c(min(x),0,max(x)),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
		axis(side=2,cex.axis=8,lwd=10,tcl=-3,padj=-1)
		mtext("Distance from feature (bp)",side=1,cex=5,outer=T,line=2)
		dev.off()

	} , mc.cores=cores, mc.preschedule=F)
	#},mc.cores=cores)
	
	# make montage
	#if(nummats < 10){
	#	system(paste("montage -geometry +2+2 -tile ",nummats,"x1 ",paste(hmnames,collapse=" ", sep=" ")," ",pwd,"/",dname,"/heatmap_montage.png",sep=""))
	#}
	
	# ##############################################################################
	# TEMPORARY IMPLEMENTATION OF CROSS-MATRIX GROUP AVERAGES
	# if(is.null(crossmat)==FALSE){
	# 	x<-((1:matcols)-(matcols/2))*winsizes[1]
	# 	pdf(file=paste(pwd,"/",dname,"/","group-averages",".pdf",sep=""))
	# 	for(i in 1:(numgroups[1])){
	# 		cat("i=",i,"\n")
	# 		groupscores=lapply(lapply(grouplist,"[",i),unlist)
	# 		allscores=unlist(groupscores)
	# 		print(allscores)
	# 		#plot(0,type="n",xlim=c(min(x),max(x)),ylim=c(min(allscores),max(allscores)),xlab="Distance from TSS",ylab="Average Score",main=paste("Group",i))
	# 		plot(0,type="n",xlim=c(min(x),max(x)),ylim=c(1,2.5),xlab="Distance from TSS",ylab="Average Score",main=paste("Group",i))
	# 		for(j in 1:length(crossmat)){
	# 			cat("j=",j,"\n")
	# 			lines(x,groupscores[[j]],lwd=3,col=matcolors[j])
	# 		}
			
	# 	}
	# 	dev.off()
	# }
	# ##############################################################################
}
# ####################################
#      KENT FUNCTIONS                #
# ####################################
getbt2index		<- function( genome ){
	if(genome=="hg19"){indexfile=paste(datapath,"hg19/igenome/Sequence/Bowtie2Index/genome",sep="")}
	if(genome=="b73v2"){indexfile=paste(datapath,"b73v2/assembly/b73v2_bowtie2/b73",sep="")}
	if(genome=="dm3"){indexfile=paste(datapath,"dm3/igenome/Sequence/Bowtie2Index/genome",sep="")}
	if(genome=="kshv"){indexfile=paste(datapath,"kshv/bt2index/kshv",sep="")}
	if(genome=="repbase"){indexfile=paste(datapath,"b73v2/repbase/repbase_bowtie2index/repbasemays",sep="")}
        return(indexfile)
}
getbtindex		<- function( genome ){
	if(genome=="b73v2"){indexfile=paste(datapath,"b73v2/btindex/genome",sep="")}
	return(indexfile)
}
getbwaindex		<- function( genome ){
	if(genome=="b73v2"){indexfile=paste(datapath,"b73v2/assembly/genome.fa",sep="")}
	return(indexfile)
}
getgenomefile		<- function( genome){
	if(genome=="b73v2"){chromfile=paste(datapath,"b73v2/b73v2.chrom.sizes",sep="")}
	if(genome=="dm3"){chromfile=paste(datapath,"dm3/dm3.chrom.sizes",sep="")}
	if(genome=="hg19"){chromfile=paste(datapath,"hg19/hg19.chrom.sizes",sep="")}
	if(genome=="hg18"){chromfile=paste(datapath,"hg18/hg18.chrom.sizes",sep="")}
	if(genome=="saccer3"){chromfile=paste(datapath,"saccer3/saccer3.chrom.sizes",sep="")}
	if(genome=="mm10"){chromfile=paste(datapath,"mm10/mm10.chrom.sizes",sep="")}
    return(chromfile)
}
gettindex		<- function( genome){
	if(genome=="hg19"){indexfile=paste(datapath,"hg19/tophatindex/genes",sep="")}
        return(indexfile)
}
getgtf			<- function( genome){
	if(genome=="hg19"){indexfile=paste(datapath,"hg19/igenome/Annotation/Genes/genes.gtf",sep="")}
        return(indexfile)
}
wigToBigWig		<- function( datafiles,genome="hg19"){
	for(l in 1:length(datafiles)){
		system(command=paste("wigToBigWig -clip", datafiles[l], getgenomefile(genome), paste(basename(removeext(datafiles[l])),".bw",sep="")))
	}
}
bedGraphToBigWig		<- function( datafiles,genome="hg19"){
	for(l in 1:length(datafiles)){
		cat(basename(datafiles[l]),": converting to bigWig\n")
		outname<-paste(basename(removeext(datafiles[l])),".bw",sep="")
		system(command=paste("bedGraphToBigWig", datafiles[l], getgenomefile(genome), outname))
		return(outname)
	}
}
bed.clip			<- function( datafiles,genome="hg19"){
	library(tools)
	for(l in 1:length(datafiles)){
		cat(basename(datafiles[l]),": clipping bed\n")
		ext<-file_ext(datafiles[l])
		outname<-paste(basename(removeext(datafiles[l])),".",ext,sep="")
		system(command=paste("bedClip", datafiles[l], getgenomefile(genome), outname))
		return(outname)
	}
}
bed.removechrom		<- function( datafiles,chrom,genome="hg19"){
	library(tools)
	for(l in 1:length(datafiles)){
		cat(basename(datafiles[l]),": removing chromosome",chrom,"from",datafiles[l],"\n")
		ext<-file_ext(datafiles[l])
		outname<-paste(basename(removeext(datafiles[l])),".",ext,sep="")
		system(command=paste("bedClip", datafiles[l], getgenomefile(genome), outname))
		return(outname)
	}
}
bedToBigBed		<- function( datafiles,genome="hg19"){
	for(l in 1:length(datafiles)){
		outname<-paste(basename(removeext(datafiles[l])),".bb",sep="")
		system(command=paste("bedToBigBed", datafiles[l], getgenomefile(genome), outname))
		return(outname)
	}
}
bigWigToBedGraph		<- function( datafiles,genome="hg19"){
	for(l in 1:length(datafiles)){
		cat(basename(datafiles[l]),": converting to bedGraph\n")
		outname<-paste(basename(removeext(datafiles[l])),".bg",sep="")
		system(command=paste("bigWigToBedGraph", datafiles[l], outname))
		return(outname)
	}
}
bigBedToBed		<- function( datafiles,genome="hg19"){
	for(l in 1:length(datafiles)){
		system(command=paste("bigBedToBed", datafiles[l], paste(basename(removeext(datafiles[l])),".bed",sep="")))
	}
}
# ####################################
#      TRACK FUNCTIONS               #
# ####################################
hub.bed			<- function( trackfiles, bedcols, color=TRUE, parentname=NULL, composite=FALSE, maxitems=1000, densecoverage=NULL, colorbystrand=NULL, tracknames=removeext(trackfiles), hubloc="dvera@epsilon.bio.fsu.edu:public_html/hubs/dlv/hg19",printtrack=F){
	library(tools)
	if(length(which(file_ext(trackfiles) %in% c("bed")) > 0)){stop("cannot create track from bed files, convert to bigBed")}
	print(data.frame("file"=trackfiles,"names"=tracknames))
	if(composite==TRUE & is.null(parentname)){stop("must specify parent name if composite=TRUE")}
	#check if dense coverage is numeric and not logical
	
	if(composite){
		track<-c(
			"",
			"",
			paste("track",parentname),
			paste("shortLabel",parentname),
			paste("longLabel",parentname),
			paste("type bigBed",bedcols,"."),
			"compositeTrack on",
			"visibility hide",
			"allButtonPair on",
			"dragAndDrop on",
			if(color){paste("itemRgb on")},
			if(is.null(colorbystrand)==FALSE){paste("colorByStrand",paste(col2rgb(colorbystrand[1]),collapse=","),paste(col2rgb(colorbystrand[2]),collapse=","),sep=" ")},
			if(is.null(densecoverage)==FALSE){paste("denseCoverage",densecoverage)},
			""
			)
	} else{
		track<-""
	}
	
	for(i in 1:length(trackfiles)){
		subtrack<-c(
			paste("track ",if(composite){paste(parentname,"-",sep="")},tracknames[i],sep=""),
			paste("bigDataUrl bbi/",trackfiles[i],sep=""),
			paste("shortLabel",tracknames[i]),
			paste("longLabel",tracknames[i]),
			paste("type bigBed",bedcols,"."),
			paste("maxItems",maxitems),
			"visilbility hide",
			if(color){paste("itemRgb on")},
			if(is.null(colorbystrand)==FALSE){paste("colorByStrand",paste(col2rgb(colorbystrand[1]),collapse=","),paste(col2rgb(colorbystrand[2]),collapse=","),sep=" ")},
			if(is.null(densecoverage)==FALSE){paste("denseCoverage",densecoverage)},
			if(composite){paste("parent",parentname)},
			""
		)
		track<-append(track,subtrack)
	}
	
	
	trackfile<-data.frame("V1"=track,stringsAsFactors=FALSE)
	if(printtrack){cat(unlist(trackfile),sep="\n")} else{write.tsv(trackfile,file="tmphubdb.txt")}
	filelist<-paste(trackfiles,sep=" ",collapse=" ")
	login<-unlist(strsplit(hubloc,":"))[1]
	path<-unlist(strsplit(hubloc,":"))[2]
	cat(paste("scp ",hubloc,"/bbi/",sep=""))
	system(paste("scp ",filelist," ",hubloc,"/bbi/",sep=""))
	system(paste("cat tmphubdb.txt | ssh ",login," 'cat >> ",path,"/hubDb.txt'",sep=""))
}

hub.bam			<- function( trackfiles, minAliQual = 0 , indelDoubleInsert = TRUE , indelQueryInsert = TRUE , bamColorMode = "gray" , bamGrayMode = "aliQual" , paired = TRUE , parentname=NULL, composite=FALSE, maxitems=1000, densecoverage=NULL, colorbystrand=NULL, tracknames=removeext(trackfiles), hubloc="dvera@epsilon.bio.fsu.edu:public_html/hubs/dlv/hg19",printtrack=F){
	library(tools)
	if(length(which(file_ext(trackfiles) %in% c("bed")) > 0)){stop("cannot create track from bed files, convert to bigBed")}
	print(data.frame("file"=trackfiles,"names"=tracknames))
	if(composite==TRUE & is.null(parentname)){stop("must specify parent name if composite=TRUE")}
	#check if dense coverage is numeric and not logical
	
	bais<-paste0(trackfiles,".bai")
	if(sum(file.exists(bais)) != length(trackfiles) ) {stop ("one or more bam indices not found!")}

	if(composite){
		track<-c(
			"",
			"",
			paste("track",parentname),
			paste("shortLabel",parentname),
			paste("longLabel",parentname),
			"type bam",
			"compositeTrack on",
			"visibility hide",
			"allButtonPair on",
			"dragAndDrop on",
			"showNames off",
			paste("minAliQual",minAliQual),
			paste("bamColorMode",bamColorMode),
			if(indelDoubleInsert){"indelDoubleInsert on"},
			if(indelQueryInsert){"indelQueryInsert on"},
			if(bamColorMode=="gray"){paste("bamGrayMode",bamGrayMode)},
			if(paired){"pairEndsByName ."},
			""
			)
	} else{
		track<-""
	}
	
	for(i in 1:length(trackfiles)){
		subtrack<-c(
			paste("track ",if(composite){paste(parentname,"-",sep="")},tracknames[i],sep=""),
			"type bam",
			if(composite){paste("parent",parentname)},
			paste("bigDataUrl bbi/",trackfiles[i],sep=""),
			paste("shortLabel",tracknames[i]),
			paste("longLabel",tracknames[i]),
			if(composite==F){paste("minAliQual",minAliQual)},
			if(composite==F){paste("bamColorMode",bamColorMode)},
			if(composite==F){"visilbility hide"},
			if(composite==F){"showNames off"},
			if(composite==F & indelDoubleInsert==T){"indelDoubleInsert on"},
			if(composite==F & indelQueryInsert==T){"indelQueryInsert"},
			if(composite==F & bamColorMode=="gray"){paste("bamGrayMode",bamGrayMode)},
			if(composite==F & paired==T){"pairEndsByName ."},
			""
		)
		track<-append(track,subtrack)
	}
	
	
	trackfile<-data.frame("V1"=track,stringsAsFactors=FALSE)
	if(printtrack){cat(unlist(trackfile),sep="\n")} else{write.tsv(trackfile,file="tmphubdb.txt")}
	filelist<-paste(trackfiles,bais,sep=" ",collapse=" ")
	login<-unlist(strsplit(hubloc,":"))[1]
	path<-unlist(strsplit(hubloc,":"))[2]
	cat(paste("scp ",hubloc,"/bbi/",sep=""))
	system(paste("scp ",filelist," ",hubloc,"/bbi/",sep=""))
	system(paste("cat tmphubdb.txt | ssh ",login," 'cat >> ",path,"/hubDb.txt'",sep=""))
}

track.wig		<- function( trackfiles,tracknames=removeext(trackfiles), plotcolors=rep("black",length(trackfiles)), user=Sys.getenv("USER"),server="epsilon.bio.fsu.edu",path="public_html/hubs/dlv",genome="hg19",range=c(0,50),printtrack=F){
	library(tools)
	if(user=="dlv04c"){user="dvera"}
	if(length(which(file_ext(trackfiles) %in% c("wig","Wig","wiggle")) > 0)){stop("cannot create track from wig files, convert to bigWig")}
	print(data.frame("file"=trackfiles,"names"=tracknames,"color"=plotcolors))
	
	track<-""
	for(i in 1:length(trackfiles)){
		c<-paste(col2rgb(plotcolors[i]),collapse=",")
		subtrack<-c(
			paste("track ",tracknames[i],sep=""),
			paste("bigDataUrl bbi/",trackfiles[i],sep=""),
			paste("shortLabel ",tracknames[i],sep=""),
			paste("longLabel ",tracknames[i],sep=""),
			paste("type bigWig ",range[1]," ",range[2],sep=""),
			paste("color ",c,sep=""),
			paste("altColor ",c,sep=""),

		""
		)
		track<-append(subtrack,track)
	}
	
	if(printtrack==TRUE){cat(unlist(outfile),sep="\n")}
	
	trackfile<-data.frame("V1"=track,stringsAsFactors=FALSE)
	write.tsv(trackfile,file="tmphubdb.txt")
	filelist<-paste(trackfiles,sep=" ",collapse=" ")
	system(paste("scp ",filelist," ",user,"@",server,":",path,"/",genome,"/bbi/",sep=""))
	system(paste("cat tmphubdb.txt | ssh ",user,"@",server," 'cat >> ",path,"/",genome,"/hubDb.txt'",sep=""))
}
hub.wig			<- function( trackfiles, range=c(0,50), parentname=NULL, multiwig=FALSE, composite=FALSE, tracknames=basename(removeext(trackfiles)), plotcolors=rainbow(length(trackfiles)), altcolors=plotcolors, hubloc="dvera@epsilon.bio.fsu.edu:public_html/hubs/dlv/hg19",printtrack=F){
	if(composite==TRUE & is.null(parentname)){stop("must specify parent name if composite=TRUE")}
	if(composite==TRUE & multiwig==TRUE){stop("cannot be both composite and multiwig")}
	print(data.frame("file"=trackfiles,"names"=tracknames,"color"=plotcolors))
	track=""
	if(multiwig){
		track<-c(
			"",
			"",
			paste("track",parentname),
			paste("shortLabel",parentname),
			paste("longLabel",parentname),
			paste("type bigWig",range[1],range[2]),
			"container multiWig",
			"graphType points",
			"configurable on",
			"visilbility hide",
			"maxHeightPixels 200:50:32",
			"aggregate transparentOverlay",
			"showSubtrackColorOnUi on",
			"autoScale on",
			"windowingFunction mean",
			"yLineOnOff on",
			"yLineMark 0",
			"smoothingWindow off",
			"alwaysZero on",
			""
			)
	}

	for(i in 1:length(trackfiles)){
		c<-paste(col2rgb(plotcolors[i]),collapse=",")
		a<-paste(col2rgb(altcolors[i]),collapse=",")
		
		subtrack<-c(
			paste("track ",if(multiwig){paste(parentname,"-",sep="")},tracknames[i],sep=""),
			paste("bigDataUrl bbi/",basename(trackfiles[i]),sep=""),
			paste("shortLabel ",tracknames[i],sep=""),
			paste("longLabel ",tracknames[i],sep=""),
			paste("type bigWig ",range[1]," ",range[2],sep=""),
			paste("color ",c,sep=""),
			paste("altColor ",a,sep=""),
			if(multiwig){paste("parent",parentname)},
			if(multiwig==F | composite==F){"graphType points"},
			if(multiwig==F | composite==F){"configurable on"},
			if(multiwig==F | composite==F){"visilbility hide"},
			if(multiwig==F | composite==F){"maxHeightPixels 200:50:32"},
			if(multiwig==F | composite==F){"autoScale on"},
			if(multiwig==F | composite==F){"windowingFunction mean"},
			if(multiwig==F | composite==F){"yLineOnOff on"},
			if(multiwig==F | composite==F){"yLineMark 0"},
			if(multiwig==F | composite==F){"smoothingWindow off"},
			if(multiwig==F | composite==F){"alwaysZero on"},
			""
		)
		track<-append(track,subtrack)
	}
	
	trackfile<-data.frame("V1"=track,stringsAsFactors=FALSE)
	if(printtrack){cat(unlist(trackfile),sep="\n")} else{write.tsv(trackfile,file="tmphubdb.txt")}
	filelist<-paste(trackfiles,sep=" ",collapse=" ")
	login<-unlist(strsplit(hubloc,":"))[1]
	path<-unlist(strsplit(hubloc,":"))[2]
	cat(paste("scp ",hubloc,"/bbi/",sep=""))
	system(paste("scp ",filelist," ",hubloc,"/bbi/",sep=""))
	system(paste("cat tmphubdb.txt | ssh ",login," 'cat >> ",path,"/hubDb.txt'",sep=""))
}
track.compositewig	<- function( trackfiles , parentname=NULL , hubloc=NULL , range=c(0,3) , plotcolors=rainbow(length(trackfiles)) , tracknames=removeext(trackfiles) ){
	
	if(is.null(parentname)){stop("must define parentname")}
	if(is.null(hubloc)){stop("must define hubloc")}
	trackfiles<-trackfiles[order(trackfiles,decreasing=T)]
	
	header<-c(
		"",
		"",
		paste("track ",parentname,"-parent",sep=""),
		"compositeTrack on",
		paste("shortLabel",parentname),
		paste("longLabel",parentname),
		"type bed 3",
		"noInherit on",
		"configurable on",
		"visilbility squish",
		paste("subGroup1 view Views ",parentname,"=",parentname,sep=""),
		"",
		paste("track",parentname),
		paste("shortLabel",parentname),
		paste("longLabel",parentname),
		paste("view ",parentname),
		paste("parent ",parentname,"-parent",sep=""),
		"graphType bar",
		"configurable on",
		"visilbility squish",
		"maxHeightPixels 200:64:32",
		"autoScale off",
		"windowingFunction mean",
		"smoothingWindow off",
		""
		)

	for(i in 1:length(trackfiles)){
		c<-paste(col2rgb(plotcolors[i]),collapse=",")
		track<-c(
			paste("track ",parentname,"-",tracknames[i],sep=""),
			paste("bigDataUrl bbi/",trackfiles[i],sep=""),
			#paste("shortLabel ",letters[i],"_",tracknames[i],sep=""),
			paste("shortLabel",tracknames[i]),
			paste("longLabel",tracknames[i]),
			paste("type bigWig",range[1],range[2]),
			paste("parent",parentname),
			paste("color",c),
			paste("altColor",c),
			paste("subGroups view=",parentname,sep=""),
			"visibility squish",
			""
		)
		header<-append(header,track)
	}
	
	outfile<-data.frame("V1"=header,stringsAsFactors=FALSE)
	write.tsv(outfile,file="tmphubdb.txt")
	filelist<-paste(trackfiles,sep=" ",collapse=" ")
	login<-unlist(strsplit(hubloc,":"))[1]
	path<-unlist(strsplit(hubloc,":"))[2]
	cat(paste("scp ",hubloc,"/bbi/",sep=""))
	system(paste("scp ",filelist," ",hubloc,"/bbi/",sep=""))
	system(paste("cat tmphubdb.txt | ssh ",login," 'cat >> ",path,"/hubDb.txt'",sep=""))
}
# ####################################
#      MISC. FUNCTIONS               #
# ####################################
isegtocoords		<- function( coords , isegs ){
	library(parallel)
	segnames<-basename(removeext(isegs))
	cat("loading coordinates\n")
	coords<-read.delim(pipe(paste("cut -f 1,2,3",coords)),stringsAsFactors=F,header=F)
	
	mclapply(1:length(isegs),function(i){
		curseg<-read.tsv(isegs[i],skip=1)
		
		mclapply(0:2, function(b){
			cat(segnames[i],": processing BC",b,"\n")
			bg<-data.frame("V1"=coords[curseg[,(b*3+1)],1],"V2"=coords[curseg[,(b*3+1)],2],"V3"=coords[curseg[,(b*3+2)],3],"V4"=curseg[,(b*3+3)])
			t<-which(complete.cases(bg$V2))
			t<-t[length(t)]
			bg<-bg[1:t,]
			outname<-paste(segnames[i],"_BC",b,".bg",sep="")
			write.tsv(bg,file=outname)
		},mc.cores=3)
	},mc.cores=floor(detectCores()/3))
}
remove.header		<- function( filename ){
	library(tools)
	ext<-file_ext(filename)
	shortname<-basename(removeext(filename))
	outname<-paste(shortname,"_rh.",ext,sep="")
	system(paste("tail -n +2",filename,">",outname))
	system(paste("mv",outname,filename))
	outname<-basename(filename)
	return(outname)
}
rpm2timing		<- function( bgfiles, reference=1, cores="max", fractionhist=TRUE , scorehist=TRUE ){
	
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	fractions=c("G1b","G2","S1","S2","S3","S4")
	
	cat("checking file count\n")
	if(floor(length(bgfiles)/6) != ceiling(length(bgfiles)/6) ){ stop("# of files not divisible by 6. missing a file?") }
	
	cat("checking file lines\n")
	if(length(unique(unlist(lapply(bgfiles,filelines)))) != 1){stop("files have different numbers of lines")}
	
	cat("finding files for each cell line\n")
	cells<-unique(remove.suffix(remove.prefix(bgfiles[grep("G2",bgfiles)],"Seq"),"G2"))
	print(cells)
	numcells<-length(cells)
	rb=rainbow(numcells)
	cellfiles<-lapply(1:numcells,function(x) bgfiles[grep(cells[x],bgfiles)] )
	print(cellfiles)
	
	outnames<-paste(cells,"_weightedscore.bg",sep="")
	finalnames<-paste(cells,"_iqr_qnorm.bg",sep="")
	#print(finalnames)
	
	cat("checking files for each cell line\n")
	numpercell=unique(unlist(lapply(cellfiles,length)))
	if(length(numpercell) != 1 | numpercell[1] != 6 ){print(cellfiles);stop("each cell line must match only 6 files")}
	checko<-lapply(1:numcells,function(x){
		lapply(1:6, function(y){
			if(grepl(fractions[y],cellfiles[[x]][y])==FALSE){
				stop(paste("cell files do not match expected fraction for",cells[x],cellfiles[[x]][y]))
			}
		} )
	})
	
	coords<-read.delim(pipe(paste("cut -f 1,2,3",cellfiles[[1]][1])),header=F)
	
	cat("calculating weighted score and IQR normalizing reference\n")
	# read in scores for each fraction for each cell line, calculated weighted score, and save file
	celldata<-mclapply(1:numcells,function(x){
		
		cd<-as.data.frame(lapply(1:6, function(y){
			as.numeric(readLines(pipe(paste("cut -f 4",cellfiles[[x]][y]))))
		} ) )
		colnames(cd)<-fractions
		return(cd)
	},mc.cores=cores )
	
	#fractiondensities<-mclapply(1:numcells,function(x){
	#	lapply(1:6,function(y){
	#		density(celldata[[x]][y])
	#	})
	#},mc.cores=cores)
	
	
	pdf(file="scoredensities.pdf")
	for(i in 1:6){
		plotmax<-quantile(celldata[[1]][,i],probs=0.95)
		plot(density(celldata[[1]][,i],from=0,to=plotmax),col=rb[1],main=fractions[i])
		for(j in 2:numcells){
			lines(density(celldata[[j]][,i],from=0,to=plotmax),col=rb[j])
		}
		legend("topright",legend=cells,col=rb,lwd=3)
	}
	
	
	celldata<-mclapply(1:numcells,function(x){
		cd<-celldata[[x]]
		wa<-(0.917*cd$G1b)+(0.750*cd$S1)+(0.583*cd$S2)+(0.417*cd$S3)+(0.250*cd$S4)+(0*cd$G2)
		wa[which(wa==0)] <- NA
		return(wa)
	},mc.cores=cores )
	
	
		plotmax<-quantile(celldata[[1]],probs=0.9,na.rm=TRUE)
		plot(density(celldata[[1]],from=0,to=plotmax),col=rb[1],main="weighted scores")
		for(j in 2:numcells){
			lines(density(celldata[[j]],from=0,to=plotmax),col=rb[j])
		}
		legend("topright",legend=cells,col=rb,lwd=3)
	dev.off()
	
	celldata<-mclapply(1:numcells,function(x){
		wa=celldata[[x]]
		if(x==reference){ wa<-( (wa-median(wa,na.rm=TRUE) ) ) * 1.59  / IQR(wa,na.rm=TRUE) }
		coords$V4=wa
		return(coords)
	},mc.cores=cores )
	
	cat("removing windows with no data\n")
	badblocks<-unique(unlist(lapply(1:numcells,function(x)(which(is.na(celldata[[x]][,4]))))))
	celldata<-mclapply(1:numcells,function(x){
		if(length(badblocks)>0){celldata[[x]]<-celldata[[x]][-badblocks,]}
		return(celldata[[x]])
	},mc.cores=cores)
	
	cat("normalizing data to",cells[reference],"and saving\n")
	celldata<-mclapply(1:numcells,function(x){
		celldata[[x]][,4][order(celldata[[x]][,4])] <- celldata[[reference]][,4][order(celldata[[reference]][,4])]
		write.tsv(celldata[[x]],file=finalnames[x])
	},mc.cores=cores)
}
rpm2timing3		<- function ( bgfiles, reference=1,cores="max"){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	
	fractions=c("G1b","G2","S1","S2","S3","S4")
	
	cat("checking file count\n")
	if(floor(length(bgfiles)/6) != ceiling(length(bgfiles)/6) ){ stop("# of files not divisible by 6. missing a file?") }
	
	#cat("checking file lines\n")
	#if(length(unique(unlist(lapply(bgfiles,filelines)))) != 1){stop("files have different numbers of lines")}
	
	cat("finding files for each cell line\n")
	cells<-unique(remove.suffix(remove.prefix(bgfiles[grep("G2",bgfiles)],"Seq"),"G2"))
	print(cells)
	numcells<-length(cells)
	rb=rainbow(numcells)
	cellfiles<-lapply(1:numcells,function(x) bgfiles[grep(cells[x],bgfiles)] )
	#print(cellfiles)
	
	ubgnames<-paste(cells,"_union.bg",sep="")
	outnames<-paste(cells,"_weightedscore.bg",sep="")
	finalnames<-paste(cells,"_iqr_qnorm.bg",sep="")
	#print(finalnames)
	
	cat("checking files for each cell line\n")
	numpercell=unique(unlist(lapply(cellfiles,length)))
	if(length(numpercell) != 1 | numpercell[1] != 6 ){print(cellfiles);stop("each cell line must match only 6 files")}
	checko<-lapply(1:numcells,function(x){
		lapply(1:6, function(y){
			if(grepl(fractions[y],cellfiles[[x]][y])==FALSE){
				stop(paste("cell files do not match expected fraction for",cells[x],cellfiles[[x]][y]))
			}
		} )
	})
	
	cat("making unionbg\n")
	#celldata<-lapply(1:numcells,function(x){
	celldata<-mclapply(1:numcells,function(x){
		cd<-read.delim(pipe(paste("bedtools unionbedg -filler NA -i",paste(cellfiles[[x]],collapse=" "))))
		colnames(cd)<-c("chrom","start","stop",fractions)
		goodrows<-which(rowSums(is.na(cd)) != 6)
		numrows<-nrow(cd)
		numgoodrows<-length(goodrows)
		cat(numrows-numgoodrows,"windows have no data in at least 1 fraction (",(numrows-numgoodrows)/numrows,"% )\n")
		cd[is.na(cd)]<-0
		wa<-(0.917*cd$G1b)+(0.750*cd$S1)+(0.583*cd$S2)+(0.417*cd$S3)+(0.250*cd$S4)+(0*cd$G2)
		cd<-data.frame("V1"=cd$chrom,"V2"=cd$start,"V3"=cd$stop,"V4"=wa,stringsAsFactors=FALSE)
		rm(wa)
		return(cd)
		print(head(cd))
	},mc.cores=cores)
	#})
	ref<-unlist(lapply(1:numcells,function(x) celldata[[x]][,4] ))
	#numrefs<-length(ref)
	
	cat("normalizing data\n")
	celldata<-mclapply(1:numcells,function(x){
		cd<-celldata[[x]]
		numpoints<-nrow(cd)
	#	if(numpoints < numrefs){
			curref<-sample(ref,numpoints)
			cd[,4][order(cd[,4])]<-curref[order(curref)]
	#	}
	#	if(numpoints > numrefs){
	#		curref<-sample(rep(ref,10),numpoints)
	#		cd[,4][order(cd[,4])]<-curref[order(curref)]
	#	}
	#	if(numpoints == numrefs){
	#		curref<-ref
	#		cd[,4][order(cd[,4])]<-curref[order(curref)]
	#	}
		cd[,4]<-( (cd[,4]-median(cd[,4],na.rm=TRUE) ) ) * 1.59  / IQR(cd[,4],na.rm=TRUE)
		write.tsv(cd,file=finalnames[x])
		return(cd)
	},mc.cores=cores)
}
scatterdens		<- function( x , y , maxscore="95%" , breaks=250 , xlims=c("1%","99%") , ylims=c("1%","99%") , white.background = F , xlabel=NULL , ylabel=NULL , cores="max" , hmcolors=c("black","blue","yellow","red") ){
	
	#check for infinite values? or just discard?
	xinfind <- which(is.infinite(x))
	yinfind <- which(is.infinite(y))
	infind <- unique(c(xinfind,yinfind))
	if(length(infind) > 0){
		x<-x[-infind]
		y<-y[-infind]
	}
	

	if(is.null(xlabel)){xlabel=deparse(substitute(x))}
	if(is.null(ylabel)){ylabel=deparse(substitute(y))}
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	rb<-colorRampPalette(hmcolors)
	#define score ranges
	if(grepl("%",xlims[1])){xlims[1]<-quantile(x,probs=as.numeric(gsub("%","",xlims[1]))/100,na.rm=TRUE)}
	if(grepl("%",xlims[2])){xlims[2]<-quantile(x,probs=as.numeric(gsub("%","",xlims[2]))/100,na.rm=TRUE)}
	if(grepl("%",ylims[1])){ylims[1]<-quantile(y,probs=as.numeric(gsub("%","",ylims[1]))/100,na.rm=TRUE)}
	if(grepl("%",ylims[2])){ylims[2]<-quantile(y,probs=as.numeric(gsub("%","",ylims[2]))/100,na.rm=TRUE)}
	
	xlims<-as.numeric(xlims)
	ylims<-as.numeric(ylims)
	
	#limit scores to defined ranges
	cat("limiting scores to range\n")
	goodrows<-which(x >= xlims[1] & x <= xlims[2] & y >= ylims[1] & y <= ylims[2])
	a<-x[goodrows]
	b<-y[goodrows]
	
	#define bin breaks
	xbreaks<-seq(xlims[1],xlims[2],(xlims[2]-xlims[1])/breaks)
	ybreaks<-seq(ylims[1],ylims[2],(ylims[2]-ylims[1])/breaks)
	
	cat("binning x-axis data\n")
	xbinind<-cut(a,xbreaks,labels=F)
	
	cat("binning y-axis data\n")
	h<-mclapply(1:breaks,function(x){
		
		hist(b[which(xbinind==x)],breaks=ybreaks,plot=F)
	},mc.cores=cores)
	
	
	binmat<-data.matrix(as.data.frame(mclapply(h,"[[",2)))
	#binmat<-binmat[nrow(binmat):1,]
	row.names(binmat)<-xbreaks[1:breaks]
	colnames(binmat)<-ybreaks[1:breaks]
	binmat<-t(binmat)
	
	if(white.background){ binmat[binmat==0]<-NA }

	cat("plotting heatmap\n")
	if(grepl("%",maxscore)){maxscore<-quantile(binmat,probs=as.numeric(gsub("%","",maxscore))/100,na.rm=TRUE)}
	if(maxscore==0){maxscore=max(binmat)}
	brks<-seq(0,maxscore,maxscore/100)
	binmat[binmat>brks[101]]<-brks[101]
	
	layout(matrix(c(4,2,1,3),nrow=2),heights=c(1,4),widths=c(4,1))
	#layout(matrix(1:2,nrow=2),heights=c(1,4))
	par(oma=c(2,2,2,2))
	par(mar=c(3,3,0,0))
	image(matrix(brks),col=rb(length(brks)-1),breaks=brks,axes=F)
	axis(side=1,at=c(0,1),labels=c(brks[1],brks[length(brks)]))
	mtext("bin counts",side=1,line=2)
	
	par(mar=c(3,3,0,0))
	image(data.matrix(binmat),breaks=brks,col=rb(length(brks)-1),xaxt='n',yaxt='n')
	axis(side=1,at=(0:4)/4,labels=xlims[1]+0:4*((xlims[2]-xlims[1])/4))
	axis(side=2,at=(0:4)/4,labels=ylims[1]+0:4*((ylims[2]-ylims[1])/4))
	mtext(xlabel,side=1,line=2)
	mtext(ylabel,side=2,line=2)
	
	dx<-density(x,from=xlims[1],to=xlims[2],na.rm=TRUE)
	dy<-density(y,from=ylims[1],to=ylims[2],na.rm=TRUE)
	par(mar=c(3,0,0,0))
	plot(dy[["y"]],dy[["x"]],type="l",yaxs='i',axes=F,lwd=4)
	par(mar=c(0,3,1,0))
	plot(dx[["x"]],dx[["y"]],type="l",xaxs='i',axes=F,lwd=4,main="")
	
}
rpm2timing2		<- function( bgfiles, reference=1, cores="max", fractionhist=TRUE , scorehist=TRUE ){
	
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	fractions=c("G1b","G2","S1","S2","S3","S4")
	
	cat("checking file count\n")
	if(floor(length(bgfiles)/6) != ceiling(length(bgfiles)/6) ){ stop("# of files not divisible by 6. missing a file?") }
	
	cat("checking file lines\n")
	if(length(unique(unlist(lapply(bgfiles,filelines)))) != 1){stop("files have different numbers of lines")}
	
	cat("finding files for each cell line\n")
	cells<-unique(remove.suffix(remove.prefix(bgfiles[grep("G2",bgfiles)],"Seq"),"G2"))
	print(cells)
	numcells<-length(cells)
	rb=rainbow(numcells)
	cellfiles<-lapply(1:numcells,function(x) bgfiles[grep(cells[x],bgfiles)] )
	print(cellfiles)
	
	outnames<-paste(cells,"_weightedscore.bg",sep="")
	finalnames<-paste(cells,"_iqr_qnorm.bg",sep="")
	#print(finalnames)
	
	cat("checking files for each cell line\n")
	numpercell=unique(unlist(lapply(cellfiles,length)))
	if(length(numpercell) != 1 | numpercell[1] != 6 ){print(cellfiles);stop("each cell line must match only 6 files")}
	checko<-lapply(1:numcells,function(x){
		lapply(1:6, function(y){
			if(grepl(fractions[y],cellfiles[[x]][y])==FALSE){
				stop(paste("cell files do not match expected fraction for",cells[x],cellfiles[[x]][y]))
			}
		} )
	})
	
	# #######################################
	
	
	coords<-read.delim(pipe(paste("cut -f 1,2,3",cellfiles[[1]][1])),header=F)
	
	cat("getting scores\n")
	# read in scores for each fraction for each cell line, calculated weighted score, and save file
	celldata<-mclapply(1:numcells,function(x){
		
		cd<-as.data.frame(lapply(1:6, function(y){
			as.numeric(readLines(pipe(paste("cut -f 4",cellfiles[[x]][y]))))
		} ) )
		colnames(cd)<-fractions
		return(cd)
	},mc.cores=cores )
	
	#fractiondensities<-mclapply(1:numcells,function(x){
	#	lapply(1:6,function(y){
	#		density(celldata[[x]][y])
	#	})
	#},mc.cores=cores)

	
	# 	pdf(file="scoredensities.pdf")
	# 	for(i in 1:6){
	# 		plotmax<-quantile(celldata[[1]][,i],probs=0.95)
	# 		plot(density(celldata[[1]][,i],from=0,to=plotmax),col=rb[1],main=fractions[i])
	# 		for(j in 2:numcells){
	# 			lines(density(celldata[[j]][,i],from=0,to=plotmax),col=rb[j])
	# 		}
	# 		legend("topright",legend=cells,col=rb,lwd=3)
	# 	}

	
	
	
	celldata<-as.data.frame(mclapply(1:numcells,function(x){
		cd<-celldata[[x]]
		wa<-(0.917*cd$G1b)+(0.750*cd$S1)+(0.583*cd$S2)+(0.417*cd$S3)+(0.250*cd$S4)+(0*cd$G2)
		wa[which(wa==0)] <- NA
		return(wa)
	},mc.cores=cores ),stringsAsFactors=FALSE)
	colnames(celldata)<-cells
	
	numwindows<-nrow(celldata)
	goodblocks<-which(complete.cases(celldata))
	cat("discarding",numwindows-length(goodblocks),"windows that have no data in at least 1 sample (",(numwindows-length(goodblocks))/numwindows,"% )\n")
	celldata<-celldata[goodblocks,]
	coords<-coords[goodblocks,]
	
	# 	plotmax<-quantile(celldata,probs=0.95,na.rm=TRUE)
	# 	plot(density(celldata[,1],from=0,to=plotmax),col=rb[1],main="weighted scores")
	# 	for(j in 2:numcells){
	# 		lines(density(celldata[,j],from=0,to=plotmax),col=rb[j])
	# 	}
	# 	legend("topright",legend=cells,col=rb,lwd=3)
	# 	dev.off()
	# 	
	reference<-rowMeans(as.data.frame(lapply(celldata,sort)))
	celldata<-as.data.frame(mclapply(1:numcells,function(x){
		celldata[,x][order(celldata[,x])]<-reference
		return(celldata[,x])
	},mc.cores=cores))
	
	celldata<-( (celldata-median(celldata,na.rm=TRUE) ) ) * 1.59  / IQR(unlist(celldata),na.rm=TRUE)
	
	celldata<-mclapply(1:numcells,function(x){
		write.tsv(cbind(coords,celldata[,x],stringsAsFactors=FALSE),file=finalnames[x])
	},mc.cores=cores)
}
qnormall<-function( seqfiles , arrayfiles, cores="max" ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	
	numseqfiles=length(seqfiles)
	numarrayfiles=length(arrayfiles)
	
	arrayheads<-lapply(1:numarrayfiles,function(x) read.delim(pipe(paste("head",arrayfiles[x]))) )
	numarraycells<-unlist(lapply(arrayheads,ncol))-2
	
	cat("pooling all data\n")
	#seqscores<-unlist(mclapply(1:numseqfiles, function(x) seqdata[[x]][,4] , mc.cores=cores))
	#arrayscores<-unlist(mclapply(1:numarrayfiles, function(x) as.vector(arraydata[[x]][,3:ncol(arraydata[[x]])]), mc.cores=cores))
	allscores<-unlist(mclapply(1:(numarrayfiles+numseqfiles), function(x){
		
		if(x<=numseqfiles){
			cat("loading",seqfiles[x],"\n")
			as.numeric(readLines(pipe(paste("cut -f 4",seqfiles[x]))))
		} else{
			a=x-numseqfiles
			cat("loading",arrayfiles[a],"\n")
			d<-unlist(read.delim(pipe(paste("cut -f",paste(3:numarraycells[a],collapse=","),arrayfiles[a]))))
			d<-as.numeric(d[2:length(d)])
			return(d)
		}
		
	}, mc.cores=cores))
	
	
	
	
	
	
	if(cores>numseqfiles){cores2=numseqfiles} else{cores2=cores}
	cat("normalizing sequencing data\n")
	mclapply(1:numseqfiles, function(x){
		seqdata<-read.tsv(seqfiles[x])
		curref<-sample(allscores,nrow(seqdata))
		seqdata[,4][order(seqdata[,4])]<-curref[order(curref)]
		outname<-paste(basename(removeext(seqfiles[x])),"_qnormToPool.bg",sep="")
		print(outname)
		write.tsv(seqdata,file=outname)
	}, mc.cores=cores2)
	
	if(cores>numarrayfiles){cores2=numarrayfiles} else{cores2=cores}
	cat("normalizing array data\n")
	mclapply(1:numarrayfiles, function(x){
		for(i in 1:numarraycells[x]){
			arraydata<-read.delim(pipe(paste("cut -f 1,2,",i+2," ",arrayfiles[x],sep="")),header=TRUE)
			curref<-sample(allscores,nrow(arraydata))
			arraydata[,3][order(arraydata[,3])]<-curref[order(curref)]
			outname<-paste(basename(removeext(arrayfiles[x])),"_",colnames(arraydata)[3],"_qnormToPool.txt",sep="")
			print(outname)
			write.tsv(arraydata,file=outname)
		}
	}, mc.cores=cores2)
}
# ####################################
#      DEATH ROW                     #
# ####################################
nolbg			<- function( bg ){
	#assumes no data within 25 bp of end of chromosome
	#assumes nol file names match column1 of bedgraph and have .txt extension
	chroms<-unique(bg$V1)
	bg$V1<-bg$V1+25
	for(i in 1:length(chroms)){
		write.tsv(bg[bg$V1==chroms[i],][,1],file="linenumbers.tmp")
		sc(paste("awk 'NR==FNR { a[$1];next } (FNR in a)' linenumbers.tmp ",chroms[i],".txt > ",chroms[i],".tmp",sep=""))
		
	}
}
mat.makenol		<- function( bedfile, regionsize=3000, start=2, stop=3, suffix="", center=FALSE, row_names=TRUE, stepsize=1,path="/lustre/maize/home/dlv04c/data/b73v2/nol/a375/"){

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
anno.dist		<- function( ){
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
vplot.make.old		<- function( fragments,bed,regionsize=1000, suffix="",ylims=c(0,200),hm.max="default",X=TRUE,savepng=TRUE){

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
vplot.make.old2		<- function( fragments,bed,regionsize=1000, suffix="",ylims=c(0,200),hm.max="default"){
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
get.file.names		<- function( paths){
	for(i in 1:length(paths)){
		namevector<-unlist(strsplit(paths[i],"/"))
		paths[i]<-namevector[length(namevector)]
	}
	paths
}
get.file.extensions	<- function( filenames){
	for(i in 1:length(filenames)){
		namevector<-unlist(strsplit(filenames[i],"\\."))
		filenames[i]<-namevector[length(namevector)]
	}
	filenames
}
# #####################################
scrapcode		<- function( ){
#if(mincovwin>0){
#				cat("finding poorly-covered features\n")
#				goodcovrows<-which(unlist(lapply(1:nrow(maskmat),function(x) length(which(maskmat[x,] > 0))))>mincovwin)
#				maskmat<-maskmat[goodcovrows,]
#				goodcovgenefile<-paste(basename(removeext(beds[i])),"_",basename(removeext(maskbed)),"_gc",mincovwin,"-",numwindows,".bed",sep="")
#				cat("saving well-covered regions in",goodcovgenefile,"\n")
#				write.tsv(bed[goodcovrows,],file=goodcovgenefile)
#
}
repliseqwindowsize	<- function( early , late , scalar="auto" , genome="hg19" , stepsize=1000 , windowsizes = 1000 * 2 ^(0:10) , cores="max"){
	
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	numwins<-length(windowsizes)
	rb<-rainbow(numwins)
	
	cat("windowing genome\n")
	windowfiles<-unlist(mclapply(1:numwins,function(x) bed.makewindows("hg19",genome=TRUE,windowsize=windowsizes[x],stepsize=stepsize,outname=paste(early,windowsizes[x],"windows",sep="")) , mc.cores=cores ))
	
	cat("counting windows\n")
	windowcounts<-unlist(mclapply(windowfiles,filelines,mc.cores=cores))
	print(windowcounts)
	
	# SORT BEDS!
	
# 	cat("calculating genome coverage\n")
# 	gcovs<-unlist(mclapply(c(early,late),bed.genomecov , genome=genome , scalar=scalar , makebigwig=FALSE , mc.cores=2 ))
# 	
# 	cat("calculating windowed coverage for early\n")
# 	ecovs<-unlist(mclapply(1:numwins,function(x) bg.window(gcovs[1], windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , operation="mean") , mc.cores=cores ))
# 	lcovs<-unlist(mclapply(1:numwins,function(x) bg.window(gcovs[2], windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , operation="mean") , mc.cores=cores ))
# 	
# 	cat("loading data\n")
# 	covs<-mclapply(ecovs,read.tsv,mc.cores=cores)
# 	
# 	covs<-mclapply(1:numwins,function(x) { 
# 		covs[[x]]$V5<-as.numeric(readLines(pipe(paste("cut -f 4",lcovs[x]))))
# 		covs[[x]]$V6<-log2(covs[[x]][,4]/covs[[x]][,5])
# 		covs[[x]]$V6[is.infinite(covs[[x]][,6])]<-NA
# 		covs[[x]]
# 	} , mc.cores=cores )
	
	cat("calculating windowed coverage for early\n")
	ecovs<-unlist(mclapply(1:numwins,function(x) bed.windowcov(early, windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , scalar=scalar) , mc.cores=cores ))
	
	cat("calculating windowed coverage for late\n")
	lcovs<-unlist(mclapply(1:numwins,function(x) bed.windowcov(late, windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , scalar=scalar) , mc.cores=cores ))
	
	ecovs<-gsub("bw","bg",ecovs)
	lcovs<-gsub("bw","bg",lcovs)
	
	cat("loading coverages\n")
	cat("loading data\n")
	covs<-mclapply(ecovs,read.tsv,mc.cores=cores)
	
	covs<-mclapply(1:numwins,function(x) { 
		covs[[x]]$V5<-as.numeric(readLines(pipe(paste("cut -f 4",lcovs[x]))))
		covs[[x]]$V6<-log2(covs[[x]][,4]/covs[[x]][,5])
		covs[[x]]$V6[is.infinite(covs[[x]][,6])]<-NA
		covs[[x]]
	} , mc.cores=cores )
	
	
	cat("counting and plotting zeroes\n")
	earlyzeroes<-(mclapply(1:numwins,function(x) which(covs[[x]][,4] == 0 ), mc.cores=cores ))
	numearlyzeroes<-unlist(lapply(earlyzeroes,length))
	
	latezeroes<-(mclapply(1:numwins,function(x) which(covs[[x]][,5] == 0 ), mc.cores=cores ))
	numlatezeroes<-unlist(lapply(latezeroes,length))
	
	bothzeroes<-(mclapply(1:numwins,function(x) which(covs[[x]][,4] == 0 & covs[[x]][,5] == 0) , mc.cores=cores ))
	numbothzeroes<-unlist(lapply(bothzeroes,length))
	
	pdf(file=paste(basename(removeext(early)),"_noscores-vs-windowsize.pdf",sep=""))
	plot(0,type="n",xlab="windowsize (bp)",ylab="% windows with no reads",xlim=c(0,128000) , ylim=c(0,50) )
	abline(v=windowsizes,col="grey70")
	points(windowsizes,100*numearlyzeroes/windowcounts,col="blue",lwd=3)
	lines(windowsizes,100*numearlyzeroes/windowcounts,col="blue",lwd=3)
	lines(windowsizes,100*numlatezeroes/windowcounts,col="red",lwd=3)
	points(windowsizes,100*numlatezeroes/windowcounts,col="red",lwd=3)
	lines(windowsizes,100*numbothzeroes/windowcounts,lwd=3)
	points(windowsizes,100*numbothzeroes/windowcounts,lwd=3)
	legend("topright",legend=c("early","late","both"),col=c("blue","red","black"),lwd=3)
	
	
	cat("calculating score distributions and plotting\n")
	equants<-unlist(mclapply(1:numwins,function(x) quantile(coverages[[x]][,4],probs=0.95,na.rm=TRUE) , mc.cores=cores ))
	equants<-1000*equants/windowsizes
	
	lquants<-unlist(mclapply(1:numwins,function(x) quantile(coverages[[x]][,5],probs=0.95,na.rm=TRUE) , mc.cores=cores ))
	lquants<-1000*lquants/windowsizes
	
	maxscore<-max(c(equants,lquants))
	
	earlydensities<-mclapply(1:numwins,function(x) density(1000*coverages[[x]][,4]/windowsizes[x],na.rm=TRUE,from=0,to=maxscore) , mc.cores=cores)
	latedensities<-mclapply(1:numwins,function(x) density(1000*coverages[[x]][,5]/windowsizes[x],na.rm=TRUE,from=0,to=maxscore) , mc.cores=cores)
	ratiodensities<-mclapply(1:numwins,function(x) density(coverages[[x]][,6],na.rm=TRUE,from=-8,to=8) , mc.cores=cores)
	
	
	plot(0,type="n",xlab="RPM",ylab="frequency of windows (kernel density estimate)",main="read density histograms for early fraction",xlim=c(0,quantile(unlist(lapply(earlydensities,"[","x")),probs=0.9)) , ylim=c(0,max(unlist(lapply(earlydensities,"[","y")))) )
	for(i in 1:numwins){
		lines(earlydensities[[i]][["x"]],earlydensities[[i]][["y"]],lwd=3,col=rb[i])
	}
	legend("topright",legend=paste(windowsizes,"bp windows") , col=rb , lwd=3, lty=1)
	
	
	plot(0,type="n",xlab="RPM",ylab="frequency of windows (kernel density estimate)",main="read density histograms for late fraction",xlim=c(0,quantile(unlist(lapply(earlydensities,"[","x")),probs=0.9)) , ylim=c(0,max(unlist(lapply(earlydensities,"[","y")))) )
	for(i in 1:numwins){
		lines(latedensities[[i]][["x"]],latedensities[[i]][["y"]],lwd=3,col=rb[i])
	}
	legend("topright",legend=paste(windowsizes,"bp windows") , col=rb , lwd=3, lty=1)
	
	
	plot(0,type="n",xlab="RT score , log2(early/late)",ylab="frequency of windows (kernel density estimate)",main="RT score histograms",xlim=c(-8,8) , ylim=c(0,max(unlist(lapply(ratiodensities,"[","y")))) )
	for(i in 1:numwins){
		lines(ratiodensities[[i]][["x"]],ratiodensities[[i]][["y"]],lwd=3,col=rb[i])
	}
	legend("topright",legend=paste(windowsizes,"bp windows") , col=rb , lwd=3, lty=1 , cex=0.75)
	
	latewherenoearly<-mclapply(1:numwins,function(x){
		coverages[[x]][,5][which(coverages[[x]][,4]==0)]
	},mc.cores=cores)
	
	earlywherenolate<-mclapply(1:numwins,function(x){
		coverages[[x]][,4][which(coverages[[x]][,5]==0)]
	},mc.cores=cores)
	dev.off()
	
	
}
	


