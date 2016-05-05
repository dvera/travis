cuffdiff2kegg <- function( gene_exp.diff , organism="hsa" , pathways="all" , limits=c(-1,1) , entrezOrganism="Hs" ){

  library(pathview)
  library(KEGGREST)
  library(gage)
  library(gageData)
  data(go.sets.hs)
  data(go.subs.hs)
  data(kegg.gs)

  de=read.tsv(gene_exp.diff,header=T)
  de$comparison=paste0(de$sample_1,"_VS_",de$sample_2)
  de$entrez=pathview::id2eg(de$gene, category ="symbol",org=entrezOrganism)[,2]
  de[which(de[,10]>10),10]<- (-10)
  de[which(de[,10]<(-10)),10]<- (10)
  de$call=0
  de$call[which(de[,14]=="yes" & de[,10]<0)] <- (-1)
  de$call[which(de[,14]=="yes" & de[,10]>0)] <- 1
  del=split(de,de$comparison)
  numcomps<-length(del)
  # save lists of differentially expressed genes
  for(i in 1:numcomps){
    unq<-as.numeric(rownames(unique(data.frame(del[[i]]$entrez))))
    del[[i]]<-del[[i]][unq,]
    row.names(del[[i]])<-del[[i]]$entrez
    write.tsv(del[[i]][,3],file=paste0(del[[i]]$comparison[1],"_allGenes.txt"))
    write.tsv(del[[i]][which(del[[i]][,14]=="yes"),3],file=paste0(del[[i]]$comparison[1],"_regulatedGenes.txt"))
    write.tsv(del[[i]][which(del[[i]][,14]=="yes" & del[[i]][10]>0),3],file=paste0(del[[i]]$comparison[1],"_upregulatedGenes.txt"))
    write.tsv(del[[i]][which(del[[i]][,14]=="yes" & del[[i]][10]<0),3],file=paste0(del[[i]]$comparison[1],"_downregulatedGenes.txt"))
  }



  # fetch target pathways
  if(length(pathways)==1 & pathways[1]=="all"){
    dbnames=keggList("pathway",organism)
    dbid=unlist(lapply(strsplit(names(dbnames),organism),"[",2))
  } else{
    if(any(grepl(organism,pathways))){
      # dbnames=unlist(lapply(keggGet(as.character(pathways)),"[",2))
      dbnames=pathways
      dbid=gsub(organism,"",pathways)
    } else{
      #dbnames=unlist(lapply(keggGet(paste0(organism,pathways)),"[",2))
      dbnames=pathways
      dbid=pathways
    }
  }


  numdbs=length(dbnames)

  # format pathway names for output files
  dbshortnames=unlist(lapply(strsplit(dbnames," - "),"[",1))
  dbshortnames=gsub(" ","",  dbshortnames)
  dbshortnames=gsub(")","",  dbshortnames)
  dbshortnames=gsub("\\(","",dbshortnames)
  dbshortnames=gsub("/","",  dbshortnames)
  dbshortnames=gsub("'","",  dbshortnames)

  for(i in 1:numcomps){

    # make a data matrix of values to preserve rownames in subsets
    vals=data.matrix(del[[i]][,c(8,9,10,13,17)])
    row.names(vals)<-del[[i]]$entrez
    # pathway enrichment analysis
    cnts.kegg.p<-gage(log2(8+vals[,1:2]),gsets=kegg.gs,ref=1,samp=2,compare="unpaired")
    o <- lapply(lapply(cnts.kegg.p,row.names),order)
    cnts.kegg.p[[1]]<-cnts.kegg.p[[1]][o[[1]],]
    cnts.kegg.p[[2]]<-cnts.kegg.p[[2]][o[[2]],]
    cnts.kegg.p[[3]]<-cnts.kegg.p[[3]][o[[3]],]

    colnames(cnts.kegg.p[[1]])<-paste0( "greater_", colnames(cnts.kegg.p[[1]]))
    colnames(cnts.kegg.p[[2]])<-paste0( "less_",    colnames(cnts.kegg.p[[2]]))
    colnames(cnts.kegg.p[[3]])<-paste0( "stats_",   colnames(cnts.kegg.p[[3]]))

    rn <- lapply(cnts.kegg.p,row.names)
    if(!identical(rn[[1]],rn[[2]]) | !identical(rn[[1]],rn[[3]])){stop("differing number of gene sets")}
    go=do.call(cbind,cnts.kegg.p)
    #write.tsv(go,file=paste0(dbshortnames[j],"_",names(del)[i]) )

    #mclapply(1:numdbs,function(j){
    for(j in 1:numdbs){
      cat(names(del)[i],": ",dbnames[j],"\n")
      pathview( gene.data=vals[,3],          pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_log2ratio_pathview"              ), sign.pos="bottomleft", kegg.native=FALSE, limit=list(cpd=limits,gene=limits) )
      pathview( gene.data=vals[,3],          pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_log2ratio_keggNative"            ), sign.pos="bottomleft", kegg.native=TRUE,  limit=list(cpd=limits,gene=limits) )
      pathview( gene.data=log10(1+vals[,1]), pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_sample1expression_pathview"      ), sign.pos="bottomleft", kegg.native=FALSE, limit=list(cpd=c(0,5),gene=c(0,5) ) )
      pathview( gene.data=log10(1+vals[,1]), pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_sample1expression_keggNative"    ), sign.pos="bottomleft", kegg.native=TRUE,  limit=list(cpd=c(0,5),gene=c(0,5) ) )
      pathview( gene.data=log10(1+vals[,2]), pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_sample2expression_pathview"      ), sign.pos="bottomleft", kegg.native=FALSE, limit=list(cpd=c(0,5),gene=c(0,5) ) )
      pathview( gene.data=log10(1+vals[,2]), pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_sample2expression_keggNative"    ), sign.pos="bottomleft", kegg.native=TRUE,  limit=list(cpd=c(0,5),gene=c(0,5) ) )
      pathview( gene.data=vals[,5],          pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_upOrDown_pathview"               ), sign.pos="bottomleft", kegg.native=FALSE, limit=list(cpd=c(-1,1),gene=c(-1,1)) )
      pathview( gene.data=vals[,5],          pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",names(del)[i], "_upOrDown_keggNative"             ), sign.pos="bottomleft", kegg.native=TRUE,  limit=list(cpd=c(-1,1),gene=c(-1,1)) )
    }
    #},mc.cores=10)
  }

  del[[1]]$color="black"
  del[[1]]$color[  which(del[[1]]$q_value < 0.1 & del[[1]]$log2.fold_change. > 0) ] = "darkgreen"
  del[[1]]$color[  which(del[[1]]$q_value < 0.1 & del[[1]]$log2.fold_change. < 0) ] = "darkred"
  del[[1]]$color[  which(del[[1]]$q_value < 0.05 & del[[1]]$log2.fold_change. > 0) ] = "green"
  del[[1]]$color[  which(del[[1]]$q_value < 0.05 & del[[1]]$log2.fold_change. < 0) ] = "red"
  gids=as.numeric(remove.prefix(read.table("http://rest.kegg.jp/link/genes/mmu00010",header=F,stringsAsFactors=FALSE)[,2],":"))
  subres = del[[1]][which(as.numeric(row.names(del[[1]])) %in% gids),]
  subres = subres[order(subres$log2.fold_change.,decreasing=TRUE),]
  barplot( subres$log2.fold_change. , names=subres$gene_id , las=3 , col=subres$color )
  #lines( 1:nrow(subres) , -log10(subres$q_value) )
  points((1:nrow(subres)*1.2)-0.5 , -log10(subres$q_value) ,col="orange",pch=20)



}
