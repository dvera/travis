vplot.make <-
function( fragmentfiles,features, regionsize=1000, windowsize=10, prunefeaturesto="", featurecenter=TRUE, strand=FALSE, start=2, stop=3, ylims=c(0,200), narrowpeak=FALSE){
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
