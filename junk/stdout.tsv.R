stdout.tsv <-
function( cmd, header=FALSE, ... ){
	read.delim(pipe(cmd),header=header,stringsAsFactors=FALSE, ... )
}
