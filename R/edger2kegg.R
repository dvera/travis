# param edgerfiles - vector of edger filepaths 
# param pathways - if not 'all', character vector of pathway ids from http://www.genome.jp/kegg/pathway.html
# param sigonly - if TRUE, only color significantly differentially expressed items, else color all 
edger2kegg <- function( edgerfiles , organism="hsa" , pathways="all" , limits=c(-1,1) , entrezOrganism="Hs", sigonly=FALSE, pval=FALSE, threads=getOption("threads",1L) ){

  library(pathview)
  library(KEGGREST)
  library(gage)
  library(gageData)
  data(go.sets.hs)
  data(go.subs.hs)
  data(kegg.gs)

  de=tsvRead(edgerfiles,col_names=T)

  # each edger file is a separate comparision, loop through each file 
  #dump <- mclapply(1:length(edgerfile),function(i) {
  for(i in 1:length(edgerfiles)) {

    # convert gene symbols to entrez IDs
    de[[i]]$entrez=pathview::id2eg(de[[i]]$gene, category ="symbol",org=entrezOrganism)[,2]

    # ???
    de[[i]][which(de[[i]]$logFC>10),10]<- (10)    # set the logFC column to 10 for each row where logFC is greater than 10
    de[[i]][which(de[[i]]$logFC<(-10)),10]<- (-10)  # set the logFC column to -10 for each row where logFC is less than -10

    # making binary calls for up/down regulated
    de[[i]]$call=0
    if(!pval) {
      de[[i]]$call[which(de[[i]]$QValue<=0.05 & de[[i]]$logFC<0)] <- (-1)  # if row is signgicant and log2(fold_change)<0, set call to -1
      de[[i]]$call[which(de[[i]]$QValue<=0.05 & de[[i]]$logFC>0)] <- 1     # if row is signgicant and log2(fold_change)>0, set call to 1
    } else {
      de[[i]]$call[which(de[[i]]$PValue<=0.05 & de[[i]]$logFC<0)] <- (-1)  # if row is signgicant and log2(fold_change)<0, set call to -1
      de[[i]]$call[which(de[[i]]$PValue<=0.05 & de[[i]]$logFC>0)] <- 1     # if row is signgicant and log2(fold_change)>0, set call to 1
    }
    unq<-as.numeric(rownames(unique(data.frame(de[[i]]$entrez)))) # get indices for rows with unique entrez ids
    de[[i]]<-de[[i]][unq,]  # remove duplicate entrez rows from de
  #### remove rows (just one i think) where entrez id is NA
    de[[i]]=de[[i]][-(which(is.na(de[[i]]$entrez))),]
  ####
    row.names(de[[i]])<-de[[i]]$entrez

    # if pathways is "all"
    if(length(pathways)==1 & pathways[1]=="all"){
      # get a list of all pathways for a species
      dbnames=keggList("pathway",organism)
      # get ID for each pathway
      dbid=unlist(lapply(strsplit(names(dbnames),organism),"[",2))
    } else{
      # if your pathway ID is species specific
      if(any(grepl(organism,pathways))){
        # dbnames=unlist(lapply(keggGet(as.character(pathways)),"[",2))
        dbnames=pathways
        # prefix pathway IDs with species code
        dbid=gsub(organism,"",pathways)
      } else{
        #dbnames=unlist(lapply(keggGet(paste0(organism,pathways)),"[",2))
        dbnames=pathways
        dbid=pathways
      }
    }

    # count number of pathways to draw
    numdbs=length(dbnames)

    # format pathway names for output files
    dbshortnames=unlist(lapply(strsplit(dbnames," - "),"[",1))
    dbshortnames=gsub(" ","",  dbshortnames)
    dbshortnames=gsub(")","",  dbshortnames)
    dbshortnames=gsub("\\(","",dbshortnames)
    dbshortnames=gsub("/","",  dbshortnames)
    dbshortnames=gsub("'","",  dbshortnames)

    # select columns of expression values
    vals=data.matrix(de[[i]][,c("logFC","logCPM","PValue","QValue","call")])
    row.names(vals)<-de[[i]]$entrez
    # if only want to color significant rows
    if(sigonly) {
      if(!pval)
        vals[which(vals[,"QValue"]>0.05),1] <- 0   # set logFC to 0 for all insignificant items (based on QValue) 
      else
        vals[which(vals[,"PValue"]>0.05),1] <- 0   # set logFC to 0 for all insignificant items (based on PValue)
    }

    dump <- mclapply(1:numdbs, function(j) {
    #for(j in 1:numdbs){
      cat(removeext(edgerfiles[i]),": ",dbnames[j],"\n")
      res1 = tryCatch({
        pathview( gene.data=vals[,1],          pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",removeext(edgerfiles[i]), "_log2ratio_pathview"              ), sign.pos="bottomleft", kegg.native=FALSE, limit=list(cpd=limits,gene=limits),low =list(gene = "red", cpd = "yellow") , mid = list(gene = "gray", cpd= "gray"), high =list(gene = "green", cpd = "blue") )
      },warning = function(w) {
          cat("\tWarning generated for pathview() call #1!\n")
      }, error = function(e) {
          cat("\tError generated for pathview() call #1!\n")
      }, finally = {
          cat("\tpathview() call #1 done.\n")
      })
      res2 = tryCatch({
        pathview( gene.data=vals[,1],          pathway.id=as.character(dbid[j]),species=organism,out.suffix=paste0(dbshortnames[j],"_",removeext(edgerfiles[i]), "_log2ratio_keggNative"            ), sign.pos="bottomleft", kegg.native=TRUE,  limit=list(cpd=limits,gene=limits),low =list(gene = "red", cpd = "yellow") , mid = list(gene = "gray", cpd= "gray"), high =list(gene = "green", cpd = "blue") )
      },warning = function(w) {
          cat("\tWarning generated for pathview() call #2!\n")
      }, error = function(e) {
          cat("\tError generated for pathview() call #2!\n")
      }, finally = {
          cat("\tpathview() call #2 done.\n")
      })
      
    }, mc.cores=threads,mc.preschedule=F) #} # \for each db
        
  } # \for each file in edgerfiles # }, mc.cores=threads,mc.preschedule=F)) 

} # \edger2kegg()

