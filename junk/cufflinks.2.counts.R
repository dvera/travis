cufflinks.2.counts <- function( genes.fpkm_tracking , samplenames=NULL){
  library(readr)
  numfiles=length(genes.fpkm_tracking)
  fl=mclapply(genes.fpkm_tracking,read_tsv,col_names=T)
  fl=mclapply(1:numfiles,function(x) fl[[x]][order(fl[[x]][,1]),] )
  e=data.matrix(as.data.frame(lapply(fl,"[",10)))
  colnames(e)=remove.suffix(f,"_")
  row.names(e)=fl[[1]][,1]

  el=log10(1+e)

}
