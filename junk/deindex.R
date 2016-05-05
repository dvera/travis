bcl2fastq <- function (cores="max"){
  cmdString <- paste(

    "bcl2fastq --create-fastq-for-index-reads",
    "-r",cores,
    "-p",cores,
    "-w",cores,
    "-d",cores,
  )

  print(cmdString)
  system(cmdString)

}

samplesheet.2.deindexsheet <- function( SampleSheet.csv, prefix ){

  csv <- read.csv(SampleSheet.csv, stringsAsFactors=F )
  csv <- csv[-(1:20),]
  i7 <- csv[,6]
  i5 <- csv[,8]
  i7u <- unique(i7)
  i5u <- unique(i5)
  i7m <- match(i7,i7u)
  i5m <- match(i5,i5u)

  idf <- data.frame(i7m,i5m,csv[,1],csv[,9],prefix)

  idf[,3]<-gsub("/","-",idf[,3])
  idf[,3]<-gsub(" ","",idf[,3])
  idf[,4]<-gsub("/","-",idf[,4])
  idf[,4]<-gsub(" ","",idf[,4])

  write.tsv(idf, file=paste0(prefix,"_sampleSheet.txt"))
  write.tsv(data.frame(i7u),file=paste0(prefix,"_i7barcodes.txt"))
  write.tsv(data.frame(i5u),file=paste0(prefix,"_i5barcodes.txt"))

}

deindex <- function( sampleSheet, i7barcodes, i5barcodes, fastqFiles, misMatch=1, cores="max" ){

  library(gtools)
  library(parallel)
  if(cores=="max"){cores=detectCores()-1}
  if(cores>length(fastqFiles)){ cores=length(fastqFiles)}

  if(any(grepl("gz",file_ext(fastqFiles)))){
    cat("gunzipping files\n")
    system(paste("echo \"",paste(fastqFiles,collapse=" "), "\" | tr ' ' '\\n' | xargs -n 1 -P",cores,"gunzip"))
    fastqFiles <- removeext(fastqFiles)
  }

  if(any(grepl("_I1_", fastqFiles))){
    cat("renaming index 1\n")
    renamefiles2(fastqFiles[grep("_I1_",fastqFiles)],"_I1_","_R3_")
    fastqFiles=gsub("_I1_","_R3_",fastqFiles)
  }

  if(any(grepl("_I2_", fastqFiles))){
    cat("renaming index 2\n")
    renamefiles2(fastqFiles[grep("_I2_",fastqFiles)],"_I2_","_R4_")
    fastqFiles=gsub("_I2_","_R4_",fastqFiles)
  }

  r1<-fastqFiles[grep("_R1_",fastqFiles)]
  r2<-fastqFiles[grep("_R2_",fastqFiles)]
  r3<-fastqFiles[grep("_R3_",fastqFiles)]
  r4<-fastqFiles[grep("_R4_",fastqFiles)]


  if(length(fastqFiles)>4){
    cat("combining files\n")
    fqs <- c("combined_R1.fastq","combined_R2.fastq","combined_R3.fastq","combined_R4.fastq")
    system(paste("cat",paste(r1,collapse=" "), "> combined_R1.fastq"))
    system(paste("cat",paste(r2,collapse=" "), "> combined_R2.fastq"))
    system(paste("cat",paste(r3,collapse=" "), "> combined_R3.fastq"))
    system(paste("cat",paste(r4,collapse=" "), "> combined_R4.fastq"))
  } else{
    fqs <- fastqFiles
  }

  cmdString <- paste(
    "deindexer -f \"rrbb\"",
    "-m", misMatch,
    "-c", sampleSheet,
    "-b", i7barcodes,
    "-b", i5barcodes,
    fqs[1],fqs[2],fqs[3],fqs[4]
  )

  print(cmdString)
  system(cmdString)

}
