bg.transform <-
function( bgfiles,operation="unlog2",cores="max"){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	numbgs<-length(bgfiles)
	bgnames<-basename(removeext(bgfiles))
	outnames<-paste(bgnames,"_",operation,".bg",sep="")
	bgs<-mclapply(bgfiles,read.tsv,mc.cores=cores)
	bgs<-mclapply(1:numbgs,function(x){
		if(operation=="unlog2"){
			bgs[[x]][,4]=2^bgs[[x]][,4]
		}
		if(operation=="log2"){
			bgs[[x]][,4]=log2(bgs[[x]][,4])
		}
		if(operation=="sqrt"){
			bgs[[x]][,4]=sqrt(bgs[[x]][,4])
		}
		if(operation=="log10" | operation=="log"){
			bgs[[x]][,4]=log(bgs[[x]][,4])
		}
		if(operation=="unlog10"){
			bgs[[x]][,4]=10^bgs[[x]][,4]
		}
		write.tsv(bgs[[x]],file=outnames[x])
	},mc.cores=cores,mc.preschedule=FALSE)
	return(outnames)
}
