shead <-
function( filename, n=10 ){
	system(paste("head -n",n,filename))
}
