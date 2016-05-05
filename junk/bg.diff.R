bg.diff <-
function( bglist1, bglist2, operation="log2", pattern="", replacement="", cores="max" ){
	if(length(bglist1) != length(bglist2)){stop("number of files in bg lists don't match")}
	numbgs<-length(bglist1)
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	if(pattern=="" | replacement==""){stop("YOU MUST SPECIFY 'pattern' AND 'replacement' FOR OUTPUT FILE NAMES")}
	outbgnames<-gsub(pattern,replacement,basename((bglist1)))
	nametab<-data.frame("BG1"=bglist1,"BG2"=bglist2,"OUTPUT"=outbgnames,stringsAsFactors=FALSE)
	print(nametab)
	if(length(which(bglist1 %in% outbgnames)==TRUE)>0){stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}
	cat("loading files\n")
	bg1<-mclapply(bglist1,read.tsv,mc.cores=cores)
	bg2<-mclapply(bglist2,read.tsv,mc.cores=cores)
	bg1lines<-unlist(lapply(bg1,nrow))
	bg2lines<-unlist(lapply(bg2,nrow))
	if(bg1lines != bg2lines){stop("files do not match in rows")}
	cat("calculating differences\n")
	bg3<-mclapply(1:numbgs,function(x){
		if(operation=="difference"){
			bg1[[x]][,4]<-bg1[[x]][,4]-bg2[[x]][,4]
		}
		if(operation=="sum"){
			bg1[[x]][,4]<-bg1[[x]][,4]+bg2[[x]][,4]
		}
		if(operation=="log2"){
			bg1[[x]][,4]<-log2(bg1[[x]][,4]/bg2[[x]][,4])
		}
		if(operation=="mean" | operation=="average"){
			bg1[[x]][,4]<-(bg1[[x]][,4]+bg2[[x]][,4])/2
		}
		write.tsv(bg1[[x]],file=outbgnames[x])
	},mc.cores=cores)
	return(outbgnames)
}
