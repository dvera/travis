cuffdiff.2.clusters <- function( gene_exp.diff, prefix, clusters=6, qval=0.15, colorder=NA , vgrep=NA ){

  library(pathview)
  library(KEGGREST)
  library(gage)
  library(gageData)
  data(go.sets.hs)
  data(go.subs.hs)
  data(kegg.gs)

  de=read.tsv(gene_exp.diff,header=T)
  if(!is.na(vgrep)){de<-de[-grep(vgrep,paste(de[,5],de[,6])),]}
  de$comparison=paste0(de$sample_1,"_VS_",de$sample_2)
  de$entrez=pathview::id2eg(de$gene, category ="symbol")[,2]
  de[which(de[,10]>10),10]<- (-10)
  de[which(de[,10]<(-10)),10]<- (10)
  de$call=0
  de$call[which(de[,14]=="yes" & de[,10]<0)] <- (-1)
  de$call[which(de[,14]=="yes" & de[,10]>0)] <- 1
  del=split(de,de$comparison)
  numcomps<-length(del)


  s1=unlist(lapply(lapply(del,"[",5),unique))
  s2=unlist(lapply(lapply(del,"[",6),unique))
  sa=c(s1,s2)
  sau=duplicated(sa)


  scores=cbind(data.frame(lapply(del,"[",8),lapply(del,"[",9)))
  row.names(scores)<-del[[1]][,3]
  scores=scores[,-which(sau)]
  colnames(scores)<-sa[-which(sau)]
  scores<-data.matrix(scores)

  qvals=data.frame(lapply(del,"[",13))
  sigs=which(apply(qvals,1,min)<qval)

  scores2=scores[sigs,]
  scores2=log2(scores2+1)
  numgenes=nrow(scores2)
  if(!is.na(colorder)){scores2<-scores2[,colorder]}
  k=kmeans(scores2,clusters)$cluster
  par(xpd=TRUE)
  image(t(scores2[order(k,decreasing=T),]),col=colorRampPalette(c("black","green"))(100))
  #segments(rep(0,numgenes),(1:numgenes)/numgenes,rep(0.05,numgenes),(1:numgenes/numgenes),lwd=3,col=rainbow(clusters)[k[order(k)]])
  lineloc=which(diff(k[order(k,decreasing=T)])!=0)/numgenes
  abline(h=lineloc,lwd=2,col="red")
  text(seq(0,1,by=1/(ncol(scores2)-1)),1.02,labels=colnames(scores2))
  knumloc<-(c(0,lineloc)+c(lineloc,1))/2
  text(1.25,knumloc,labels=clusters:1,cex=2)


  scores3<-as.data.frame(scores2)
  scores3$cluster<-k
  scores4<-split(scores3,scores3$cluster)

  X11()

  plot(0,type="n",xlim=c(1,1+(ncol(scores2))),ylim=c(0,ceiling(max(scores2))))
  rb=rainbow(clusters)

  for(i in 1:clusters){
    lines(1:ncol(scores2),colMeans(scores4[[i]][,1:ncol(scores2)]),type="l",col=rb[i],lwd=5)
  }
  legend("topright",legend=1:clusters,col=rb,lwd=5)

  write.tsv(scores3,file=paste0(prefix,"_clusters.tsv"))




}
