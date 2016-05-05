stail <-
function( filename, n=10 ){
	system(paste("tail -n",n,filename))
}
