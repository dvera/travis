qteller.2.matrix <-
function (cuffdiff.file , template.matrix.file ,  b73=FALSE){
	library(tools)
	matsuffix<-get.suffix(basename(template.matrix.file),"_")
	winsize=as.numeric(gsub("mat","",file_ext(template.matrix.file)))

	mat<-read.mat(template.matrix.file)
	cuf<-read.csv(cuffdiff.file,header=T)
	matgenes<-unlist(lapply(strsplit(row.names(mat),";") , "[" , 3 ))
	if(b73){ 
		matgenes<-remove.suffix(matgenes,"_T")
		matgenes<-gsub("_FGT","_FG",matgenes)
	}
	namematch<-match(matgenes,cuf[,1])
	 
	
	for(i in 6:(ncol(cuf)-1 ) )  {
		samplename<-colnames(cuf)[i]
		outname <- paste0 (samplename,"_FPKM_",matsuffix)
		m1<-matrix(as.numeric(cuf[namematch,i]),nrow=nrow(mat),ncol=ncol(mat))
		row.names(m1)<-row.names(mat)
		write.mat(m1,file=outname)
	
	}	
}
