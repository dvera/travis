bg.uniondiff <-
function( bglist1, bglist2, pattern="", replacement="", cores="max" ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	if(pattern=="" | replacement==""){stop("YOU MUST SPECIFY 'pattern' AND 'replacement' FOR OUTPUT FILE NAMES")}
	bglist1names<-basename(bglist1)
	bglist2names<-basename(bglist2)
	outbgnames<-gsub(pattern,replacement,bglist1names)
	nametab<-data.frame("BG1"=bglist1names,"BG2"=bglist2names,"OUTPUT"=outbgnames,stringsAsFactors=FALSE)
	print(nametab)
	if(length(which(bglist1names %in% outbgnames)==TRUE)>0){stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}
	cat("making unionbg\n")
	system(paste("bedtools unionbedg -i",paste(c(bglist1,bglist2),collapse=" "),"> tmpunionbg.tmp"))
	#cat(paste("bedtools unionbedg -i",paste(c(bglist1,bglist2),collapse=" "),"> tmpunionbg.tmp\n"))
	cat("calculating differences\n")
	mclapply(1:length(bglist1), function(x){
		#cat(paste("awk '{print $1,$2,$3,$",3+x,"-$",3+x+length(bglist1),"' OFS='\\t' tmpunionbg.tmp | cut -f 1,2,3,4 > ",outbgnames[x],"\n",sep=""))
		system(paste("awk '{print $1,$2,$3,$",3+x,"-$",3+x+length(bglist1),"}' OFS='\\t' tmpunionbg.tmp | cut -f 1,2,3,4 > ",outbgnames[x],sep=""))
	},mc.cores=cores)
}
