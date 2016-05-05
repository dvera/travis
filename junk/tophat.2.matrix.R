tophat.2.matrix <-
function (genes.fpkm_tracking , templateMatrix , samplename , b73=FALSE){

  library(tools)

  matsuffix<-get.suffix(basename(templateMatrix),"_")
	winsize=as.numeric(gsub("mat","",file_ext(templateMatrix)))

	mat<-read.mat(templateMatrix)
	expdata<-read.tsv(genes.fpkm_tracking,header=T)

	matgenes<-unlist(lapply(strsplit(row.names(mat),";") , "[" , 3 ))

  if(b73){
		matgenes<-remove.suffix(matgenes,"_T")
		matgenes<-gsub("_FGT","_FG",matgenes)
	}

	namematch<-match(matgenes,expdata[,1])

	m<-matrix(as.numeric(expdata[namematch,8]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(m1)<-row.names(mat)

	m1name<-paste0(samplename,"_FPKM_",matsuffix)

	write.mat(m,file=m1name)

  return(m1name)

}
