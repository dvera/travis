cuffdiff2list<-function(gene_exp.diff,pval=0.05){
  library(pathview)
  library(KEGGREST)
  library(gage)
  library(gageData)
  data(go.sets.hs)
  data(go.subs.hs)
  data(kegg.gs)

  de=read.tsv(gene_exp.diff,header=T)
  de$comparison=paste0(de$sample_1,"_VS_",de$sample_2)
  de$entrez=pathview::id2eg(de$gene, category ="symbol")[,2]
  # de[which(de[,10]>10),10]<- (-10)
  # de[which(de[,10]<(-10)),10]<- (10)
  de$call=0
  de$call[which(de[,13] < pval & de[,10]<0)] <- (-1)
  de$call[which(de[,13] < pval & de[,10]>0)] <- 1
  del=split(de,de$comparison)
  numcomps<-length(del)
  # save lists of differentially expressed genes
  for(i in 1:numcomps){
    unq<-as.numeric(rownames(unique(data.frame(del[[i]]$entrez))))
    del[[i]]<-del[[i]][unq,]
    #row.names(del[[i]])<-del[[i]]$entrez
    write.tsv(del[[i]][,3],file=paste0(del[[i]]$comparison[1],"_allGenes.txt"))
    write.tsv(del[[i]][which(del[[i]]$call != 0  ),3],file=paste0(del[[i]]$comparison[1],"_regulatedGenes.txt"))
    write.tsv(del[[i]][which(del[[i]]$call == 1  ),3],file=paste0(del[[i]]$comparison[1],"_upregulatedGenes.txt"))
    write.tsv(del[[i]][which(del[[i]]$call == -1 ),3],file=paste0(del[[i]]$comparison[1],"_downregulatedGenes.txt"))
  }


}
