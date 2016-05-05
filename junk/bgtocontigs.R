bgtocontigs <-
function(bgfile, cores="max", maxgap=1 ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	cat("loading bedgraph\n")
	bg<-read.tsv(bgfile)
	stepsize=unique(bg$V3-bg$V2)
	if(length(stepsize)>1){stop("ERROR: bg interval size variation exists")} else{cat("step size =",stepsize)}
	cat("splitting pos/neg scores\n")
	p<-bg
	n<-bg
	rm(bg)
	samplename<-basename(removeext(bgfile))
	p$V4[which(p$V4<0)]<-0
	n$V4[which(n$V4>0)]<-0

	cat("finding gaps\n")
	gaps<-p[2:nrow(p),2]-p[1:(nrow(p)-1),3]
	stops<-which( gaps > maxgap )
	starts=c(1,stops+1)
	stops=c(stops,nrow(p))
	dir.create("pos")
	dir.create("neg")
	dir.create("int")
	cat("saving data\n")
	dir.create(samplename)
	dir.create(paste0(samplename,"/pos"))
	dir.create(paste0(samplename,"/neg"))
	dir.create(paste0(samplename,"/int"))
	a<-mclapply(1:length(starts),function(x){
		write.tsv(p[starts[x]:stops[x],4],file=paste(samplename,"/pos/",paste(p[starts[x],1:2],collapse="-"),".pos",sep=""))
		write.tsv(n[starts[x]:stops[x],4],file=paste(samplename,"/neg/",paste(n[starts[x],1:2],collapse="-"),".neg",sep=""))
		write.tsv(n[starts[x]:stops[x],1:3],file=paste(samplename,"/int/",paste(n[starts[x],1:2],collapse="-"),".int",sep=""))
	},mc.cores=cores)
	write.tsv(p,file=paste(basename(removeext(bgfile)),"_pos.bg",sep=""))
	write.tsv(n,file=paste(basename(removeext(bgfile)),"_neg.bg",sep=""))
}
