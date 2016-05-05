isegtobed <-
function(posfile,negfile=NULL,intervalfiles,prefix="z3b1sens",cutoffs=c(0,1,1.5,2,3),genome="b73v2", fix=FALSE){
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
