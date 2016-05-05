isegtocoords <-
function( coords , isegs ){
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
