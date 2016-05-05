bgCorMatrix <- function( bgs , metadata=NULL , cores=2 , cormethod="pearson", threads=getOption("threads",1L), ... ){

  #library(pheatmap)
  if(!is.null(metadata)){md=read.tsv(metadata,header=TRUE,row.names=1)}

  numfiles=length(bgs)

  fl=bgRead(bgs, threads=threads, makematrix=TRUE , enforceEquality=TRUE)
  fl <- as.data.frame(fl)
  eg=expand.grid(1:numfiles,1:numfiles)

  pb <- txtProgressBar(min = 0, max = numfiles*numfiles, style = 3)

  lc=as.data.frame(
    matrix(
      unlist(mclapply(1:nrow(eg),function(x){
        val=cor(fl[[eg[x,1]]],fl[[eg[x,2]]],method=cormethod)
        setTxtProgressBar(pb, x)
        return(val)
      },mc.cores=threads,mc.preschedule=T)),
    nrow=numfiles)
  )
  cat("\n")

  #fheatmap(lcm,annotation_col=ac,cluster_rows=T,display_tree_row=T,annotation_row=ar)
  rownames(lc) <- basename(removeext(bgs))
  colnames(lc) <- rownames(lc)
  #fheatmap(lc,cluster_rows=T,cluster_cols=T,display_tree_col=T,display_tree_row=T,names_font_style="mono")
  #pheatmap(lc, annotation_row=md, annotation_col=md , fontsize_row=5, fontsize_col=5,  ...)

  # if(!is.null(metadata)){
  #   factorcols=which(unlist(lapply(md,class))=="factor")
  #   factorcolors=lapply(1:length(fcs),function(x) cols[1:length(table(md[factorcols[x]]))])
  #   for(i in 1:length(factorcols)){
  #     names(factorcolors[[i]])<-unique(as.character(md[,factorcols[i]]))
  #   }
  #   names(factorcolors)<-colnames(md)[factorcols]
  # }
  # pheatmap(
  #   lc,
  #   annotation_row=if(is.null(metadata)){NA} else{md[g1,-(9:13)]},
  #   annotation_col=if(is.null(metadata)){NA} else{md[g1,-(9:13)]} ,
  #   fontsize_row=10,
  #   fontsize_col=10,
  #   annotation_colors=if(is.null(metadata)){NA} else{factorcolors}
  # )
  return(lc)
}
