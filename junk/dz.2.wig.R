dz.2.wig <-
function( dzfile, cores="max", splitchromosomes=FALSE ){
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
