bedCat <- function (beds , outname , buffersize="10G" , cores=8 , presorted=TRUE ){

	cmdString <- paste(
		"sort -T . -S",
		buffersize,
		"--parallel",cores,
		if(presorted){"-m"},
		"-k1,1 -k2,2n",
		"-o",
		outname,
		paste(beds,collapse=" ")
	)

	print(cmdString)
	system(cmdString)
	return(outname)

}
