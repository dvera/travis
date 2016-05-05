bg.dist <-
function( bgfile, bedfiles, cutoff = 1, piefeatures=c("intergenic","promoter5","5utr","CDS","3utr","promoter3"),cores="max" ){
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
