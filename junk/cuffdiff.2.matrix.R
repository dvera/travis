cuffdiff.2.matrix <-
function (cuffdiff.file , template.matrix.file , sample1name , sample2name , b73=FALSE){
	library(tools)
	matsuffix<-get.suffix(basename(template.matrix.file),"_")
	winsize=as.numeric(gsub("mat","",file_ext(template.matrix.file)))

	mat<-read.mat(template.matrix.file)
	cuf<-read.tsv(cuffdiff.file,header=T)
	matgenes<-unlist(lapply(strsplit(row.names(mat),";") , "[" , 3 ))
	if(b73){ 
		matgenes<-remove.suffix(matgenes,"_T")
		matgenes<-gsub("_FGT","_FG",matgenes)
	}
	namematch<-match(matgenes,cuf[,1])
	 
	m1<-matrix(as.numeric(cuf[namematch,8]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(m1)<-row.names(mat)
	
	m2<-matrix(as.numeric(cuf[namematch,9]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(m2)<-row.names(mat)

	cuf$invq<-1/cuf[,13]
	gi2<-which(cuf[,9]-cuf[,8]>0)
	cuf$invq[gi2]<-cuf$invq[gi2]*-1
	invq<-matrix(as.numeric(cuf$invq[namematch]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(invq)<-row.names(mat)

	dif<-matrix(as.numeric(cuf[namematch,8]-cuf[namematch,9]),nrow=nrow(mat),ncol=ncol(mat))
	row.names(dif)<-row.names(mat)

	invqname<-paste0(sample1name,"-",sample2name,"_inverseQval_",get.suffix(basename(template.matrix.file),"_"))
	m1name<-paste0(sample1name,"_FPKM_",matsuffix)
	m2name<-paste0(sample2name,"_FPKM_",matsuffix)
	difname<-paste0(sample1name,"Minus",sample2name,"_FPKM_",matsuffix)

	write.mat(m1,file=m1name)
	write.mat(m2,file=m2name)
	write.mat(dif,file=difname)
	write.mat(invq,file=invqname)
	
}
