bedStructures <-
function( bedfile , promoter5 = c(-1000,0) , promoter3 = c(0,1000) , bedname = basename(removeext(bedfile)), chromsizes ) {

	if(missing(chromsizes)){
		chromsizes<-getOption("chromsizes",NULL)
		if(is.null(chromsizes)){stop("must define file contain chromosome sizes")}
	}
	if(!file.exists(chromsizes)){
		stop("chromsizes file does not exist")
	}

	utr5name  <- paste0(bedname,"_utr5.bed")
	utr3name  <- paste0(bedname,"_utr3.bed")
	utrname   <- paste0(bedname,"_utr.bed")
	prom5name <- paste0(bedname,"_prom5.bed")
	prom3name <- paste0(bedname,"_prom3.bed")
	cdsname   <- paste0(bedname,"_cds.bed")
	exonname  <- paste0(bedname,"_exons.bed")
	intergenicname <- paste0(bedname,"_intergenic.bed")
	genicname <- paste0(bedname,"_genic_tmp.bed")
	genicname2 <- paste0(bedname,"_genic.bed")
	intronname <- paste0(bedname,"_introns.bed")

	bed<-read.table(bedfile,header=F,stringsAsFactors=T,sep="\t")

	bed<-bed[order(bed[,1],bed[,2],bed[,3]),]

	if(ncol(bed) < 6){stop("no strand information found")}
	posbed<-bed[which(bed[,6]=="+"),]
	negbed<-bed[which(bed[,6]=="-"),]

	chroms = unique(bed[,1])
	sbed=split(bed,bed[,1])
	numchroms=length(sbed)

	exons <- do.call(rbind,lapply(1:numchroms,function(i){
		exonsizes <- lapply(strsplit(sbed[[i]][,11],","),as.numeric)
		exonstarts <- lapply(strsplit(sbed[[i]][,12],","),as.numeric)
		exonstarts2 <- unlist(lapply(1:length(exonsizes),function(x) exonstarts[[x]]+sbed[[i]][x,2] ))
		exonsizes2 <- unlist(exonsizes)
		exonstops <-  exonstarts2+exonsizes2
		exons <- data.frame(chrom=names(s)[[i]],start=exonstarts2,stop=exonstops)
		return(exons)
	}))

	tsvWrite(exons,file=exonname)
	bedSort(exonname)

	utr5<-rbind(
		data.frame(V1=posbed[,1],V2=posbed[,2],V3=posbed[,7],stringsAsFactors=F),
		data.frame(V1=negbed[,1],V2=negbed[,8],V3=negbed[,3],stringsAsFactors=F)
	)
	utr5<-utr5[which(utr5[,3]-utr5[,2] > 0),]
	tsvWrite(utr5,file=utr5name)
	bedSort(utr5name)

	utr3<-rbind(
		data.frame(V1=posbed[,1],V2=posbed[,8],V3=posbed[,3],stringsAsFactors=F),
		data.frame(V1=negbed[,1],V2=negbed[,2],V3=negbed[,7],stringsAsFactors=F)
	)
	utr3<-utr3[which(utr3[,3]-utr3[,2] > 0),]
	tsvWrite(utr3,file=utr3name)
	bedSort(utr3name)

	utrname <- bedCat(c(utr5name,utr3name),utrname)

	prom5<-rbind(
		data.frame(V1=posbed[,1],V2=posbed[,2]+promoter5[1],V3=posbed[,2]+promoter5[2],stringsAsFactors=F),
		data.frame(V1=negbed[,1],V2=negbed[,3]-promoter5[2],V3=negbed[,3]-promoter5[1],stringsAsFactors=F)
	)
	prom5[which(prom5[,2]<0),2]<-0
	prom5=prom5[which((prom5[,3]-prom5[,2])>0),]
	tsvWrite(prom5,file=prom5name)
	bedSort(prom5name)
	prom5name <- kentBedClip(prom5name,chromsizes)

	prom3<-rbind(
		data.frame(V1=posbed[,1],V2=posbed[,3]+promoter3[1],V3=posbed[,3]+promoter3[2],stringsAsFactors=F),
		data.frame(V1=negbed[,1],V2=negbed[,2]-promoter3[2],V3=negbed[,2]-promoter3[1],stringsAsFactors=F)
	)
	prom3[which(prom3[,2]<0),2]<-0
	prom3=prom3[which((prom3[,3]-prom3[,2])>0),]
	tsvWrite(prom3,file=prom3name)
	bedSort(prom3name)
	prom3name <- kentBedClip(prom3name,chromsizes)


	bedfile2name=paste0(bedfile,"_gene.bed")
	cmdString <- paste("cut -f 1,2,3",bedfile,">",bedfile2name)
	res <- cmdRun(cmdString)

	intronname <- bedtoolsSubtract(bedfile2name,exonname,intronname)
	bedSort(intronname)
	cdsBoundaries <- bedtoolsSubtract(bedfile2name,utrname)
	bedSort(cdsBoundaries)
	cdsname <- bedtoolsSubtract(cdsBoundaries,intronname,cdsname)
	bedSort(cdsname)

	genicname=bedCat(c(bedfile2name,prom5name,prom3name),genicname)
	genicname=bedtoolsMerge(genicname)
	cmdString <- paste("cut -f 1,2,3",genicname,">",genicname2)
	res <- cmdRun(cmdString)
	intergenicfile <- bedtoolsComplement(genicname2,chromsizes,intergenicname)
	unlink(cdsBoundaries)
	unlink(genicname)
	unlink(genicname2)



	beds=c(intergenicname,intronname,utr5name,utr3name,prom5name,prom3name,cdsname)
	beds=kentBedClip(beds,chromsizes)

}
