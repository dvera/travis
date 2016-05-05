downloadfile <-
function( url){
	if(Sys.info()["sysname"]=="Linux"){
		sc(paste("wget -q ",url,sep=""))
	}
	if(Sys.info()["sysname"]=="Darwin"){
		sc(paste("curl -silent -O ",url,sep=""))
	}
}
