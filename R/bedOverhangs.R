bedOverhangs <- function( bedfile, chrom, start, stop, distrange=c(-20,20), cores="max" ){

	library(parallel)
	if(cores=="max"){cores<-detectCores()-1}

	cat("loading bed file\n")
	bed<-read.delim(pipe(paste0("awk '{ if($1==\"",chrom,"\" && $2>=",start," && $3<=",stop,") print $0 }' ", bedfile)), stringsAsFactors=FALSE, header=F)

	ps<-bed[which(bed$V6=="+"),]
	ns<-bed[which(bed$V6=="-"),]


	d<-unlist(mclapply(1:nrow(ps), function(x){
		posend <- ps[x,2]
		lflank <- posend+distrange[1]
		rflank <- posend+distrange[2]
		negend<-ns[,3]
		negs <- ns[which(negend>lflank & negend<rflank),2]
		distances<-posend-negs+1
		return(distances)
	},mc.cores=cores,mc.preschedule=T))


	return(d)






}
