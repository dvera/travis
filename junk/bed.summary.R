bed.summary <- function( bedfiles , genes , querybed , genomefile , outname , slop=1000 , bednames=basename(removeext(bedfiles)) , cores="max" ){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	if(length(bedfiles) < cores){cores=length(bedfiles)}
	numbeds <- length(bedfiles)
	numgenes <- filelines(genes)
	qbed <- querybed
	qbed <- bedtoolsMerge(qbed)
	totalBp <- bed.breadth(qbed)
	counts <- unlist(mclapply(bedfiles, filelines, mc.cores=cores))
	breadth <- unlist(mclapply(bedfiles, bed.breadth, mc.cores=cores))
	percentOftotalBp <- 100*breadth/totalBp
	numPerKb <- 1000*counts/totalBp
	slopgenes <- bed.slop(genes,genomefile,slop,slop)
	neargenebeds <- unlist(mclapply(bedfiles,bedtoolsIntersect,b=slopgenes,mc.cores=cores))
	numneargenes <- unlist(mclapply(neargenebeds,filelines,mc.cores=cores))
	pctNearGenes <- numneargenes/counts
	nearbedgenes <- unlist(mclapply(1:numbeds, function(x) bedtoolsIntersect(slopgenes,bedfiles[x]),mc.cores=cores))
	numgenesnear <- unlist(mclapply(nearbedgenes,filelines,mc.cores=cores))
	pctGenesNear <- 100* numgenesnear / numgenes
	medianSize <- as.numeric(
		unlist(mclapply(bedfiles, function(x)
			system(paste("awk '{print $3-$2}'",x," | sort -k1,1n | awk 'BEGIN{i=0} {a[i]=$1 ; i++} END{print a[int(NR/2)]}'"),intern=TRUE),
			mc.cores=cores))
		)

	medianSizeNearGenes <- as.numeric(
		unlist(mclapply(neargenebeds, function(x)
			system(paste("awk '{print $3-$2}'",x," | sort -k1,1n | awk 'BEGIN{i=0} {a[i]=$1 ; i++} END{print a[int(NR/2)]}'"),intern=TRUE),
			mc.cores=cores))
		)
	
	summarytable <- data.frame(bednames, counts, totalBp, breadth, percentOftotalBp, numPerKb, pctNearGenes, pctGenesNear, medianSize, medianSizeNearGenes)

	write.table( summarytable , file=outname , sep="\t" , col.names=T, row.names=F, quote=F )

	return(outname)
}