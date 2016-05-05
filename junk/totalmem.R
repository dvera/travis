totalmem <-
function( ){
	return(as.numeric(readLines(pipe("free -b | awk '{print $2}' | head -n 2 | tail -n 1"))))
}
