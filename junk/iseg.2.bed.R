isegToBed <- function ( isegfiles , intfiles , cores = 11 ){
	intnames<-basename(removeext(intfiles))
	isegnames<-basename(removeext(isegfiles))
	
	cat("reading in interval files\n")
	int<-mclapply(intfiles,read.tsv,mc.cores=cores,mc.preschedule=F)
	
	cat("reading in iseg files\n")
	iseg<-mclapply(isegfiles,read.csv,header=T,stringsAsFactors=F,mc.cores=cores,mc.preschedule=F)
	
	beds<-mclapply(1:length(iseg), function(x){
		cat("processing",isegnames[x],"\n")
		bedrows <- lapply(1:nrow(iseg[[x]]), function(y){
			intmatch<-which(intnames==iseg[[x]][y,1])
			c( iseg[[x]][y,2] , int[[intmatch]][iseg[[x]][y,4],2] , int[[intmatch]][iseg[[x]][y,5],3] )
		})
		bedrows <- as.data.frame(do.call(rbind,bedrows),stringsAsFactors=FALSE)
		bedrows$V2<-as.numeric(bedrows$V2)
		bedrows$V3<-as.numeric(bedrows$V3)
		bedrows<-bedrows[which(bedrows$V3-bedrows$V2 > 0) , ]
		bedrows <- bedrows [ order(bedrows$V1,bedrows$V2) , ]
		write.tsv(bedrows, file = paste0 ( isegnames[x] , ".bed"))
	}, mc.cores = cores, mc.preschedule = F)

}