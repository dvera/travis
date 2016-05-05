vplot.rank <-
function( peakfiles, bgfiles, fragments, cellline="Gm12878", cores=16, np=TRUE, ... ){
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
