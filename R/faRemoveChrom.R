faRemoveChrom <- function( fastaFiles , removePatterns , suffix="_rmChr" ){

  
  outnames<-paste0(basename(removeext(fastaFiles)),suffix,".fa")
	searchString <- paste0(paste0("$1 ~ /",removePatterns),"/",collapse=" || ")
  cmdString <- paste("awk '{if($1~/^>/){if(",searchString,"){a=0}else{a=1}};if(a==1){print $0}}'",fastaFiles,">",outnames)
  
}
