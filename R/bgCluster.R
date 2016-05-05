bgCluster <- function( bedGraphs , threshold=NULL , legendnames=NULL , gradient=c("red","black","green") , breaks=100 , range=NULL  ){

  scores <- bgRead(bedGraphs, makematrix=T, enforceEquality=T, bgnames=legendnames, threads=12)
  bg <- tsvRead(bedGraphs[1])
  rownames(scores) <- paste(bg[,1],bg[,2],bg[,3],sep="-")

  if( !is.null(threshold) & is.numeric(threshold) ){
    gt <- which(apply(scores,1,max) > threshold )
    lt <- which(apply(scores,1,min) < -threshold )
    glt <- unique(c(gt,lt))
    scores <- scores[glt,]
  }

  if(is.null(range)){
    range=c(min(scores),max(scores))
  } else{
    scores[scores>range[2]] <- range[2]
    scores[scores<range[1]] <- range[1]
  }
  heatmap.2(t(data.matrix(scores)), col=colorRampPalette(gradient)(breaks), breaks=seq(range[1],range[2],length.out=breaks+1),trace='none' ,margins=c(20,20))
}
