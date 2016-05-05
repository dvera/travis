bgScatterGrid <-
function( bgfiles, bgnames=NULL, threads=getOption("threads",1L), ...  ){

	#bgnames<-basename(removeext(bgfiles))
	numbgs<-length(bgfiles)

	# noplots <- unlist(lapply(1:numbgs, function(i)  print( (i-1)*numbgs+(i:numbgs) ) ))
	# dgl <- diag(matrix(1:(numbgs^2),nrow=numbgs))

	# colramp <- colorRampPalette(c("red","black","green"), space = "rgb")
	# brks=seq(-1,1,by=2/100)
	# hcolors<-colramp(100)

	#make list of matrices
	#cat("reading in bedgraphs\n")
	if(is.null(bgnames)){ bgnames <- basename(removeext(bgfiles)) }

	bglist<-bgRead(bgfiles,threads=threads)

	#make pairwise matrix correlations
	#cat("calculating global correlations\n")
	pairs<-expand.grid(1:numbgs,1:numbgs)
	par(mfrow=c(numbgs,numbgs),oma=c(4,4,4,4),mar=c(0.2,0.2,0.2,0.2))
	bottoms <- numbgs^2:(numbgs^2-numbgs+1)
	lefts <- seq(1,numbgs^2,numbgs)
	for(i in 1:(numbgs^2 ) ) {
		scatterdens(bglist[pairs[i,1],] , bglist[pairs[i,2],],basic=T, ... , xaxislabel=if(i %in% bottoms){T} else{F}, yaxislabel=if(i %in% lefts){T} else{F}, xlabel=bgnames[pairs[i,1]], ylabel=bgnames[pairs[i,2]] )
	}
	#bgcors<-as.numeric(unlist(mclapply(1:nrow(pairs),function(x){
	#	cor(bglist[[pairs[x,1]]] , bglist[[pairs[x,2]]] , use="complete.obs", method="spearman" )
	#},mc.cores=cores)))
	#bgcors<-matrix(bgcors,nrow=numbgs,ncol=numbgs)
	#row.names(bgcors)<-bgnames
	#colnames(bgcors)<-bgnames



	#}
	#globaltablename<-paste(matnames[1],"_globalcors.tsv",sep="")
	#cat("saving global correlation table to",globaltablename,"\n")
	#write.table(cormat,file=globaltablename,sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)





}
