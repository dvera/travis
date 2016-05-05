bam.parseBySnp <- function ( erlyBams , lateBams , vcfFile , sampleNames , genomefile, include=NULL , exclude=NULL , quiet=TRUE , cores="max" ){

  library(parallel)
  if(cores=="max"){cores<-detectCores()-1}
  numbams=length(c(erlyBams))
  if(numbams<cores){cores<-numbams}

  vcf=vcfFile

  if(!all(file.exists(paste0(c(erlyBams,lateBams),".bai")))){stop("bams are not indexed")}

  erlyVarSupportSams <- paste0(basename(removeext(erlyBams)),"_supportingVariant_",basename(removeext(vcf)),".sam")
  erlyRefSupportSams <- paste0(basename(removeext(erlyBams)),"_supportingReference_",basename(removeext(vcf)),".sam")
  lateVarSupportSams <- paste0(basename(removeext(lateBams)),"_supportingVariant_",basename(removeext(vcf)),".sam")
  lateRefSupportSams <- paste0(basename(removeext(lateBams)),"_supportingReference_",basename(removeext(vcf)),".sam")

  cat("counting reads\n")
  erlyCounts<-unlist(mclapply(erlyBams,bam.count,mc.cores=cores))
  lateCounts<-unlist(mclapply(lateBams,bam.count,mc.cores=cores))
  erlyFactors<- 1000000/erlyCounts
  lateFactors<-  1000000/lateCounts

  if(!is.null(include)){
    vcf=bed.intersect(vcf,include)
  }
  if(!is.null(exclude)){
    vcf=bed.intersect(vcf,include,extraargs="-v")
  }

  cat("parsing reads\n")
  cmdString1 <- paste0(
    "vcf_read_support.py ",
    " -b ", erlyBams,
    " -i ", vcf,
    if(quiet){" 2>/dev/null "},
    #" | awk '{ ",
    " | awk '{ if( !a[$1]++ ){",
      " if($0 ~ \"\tvr\"){print $0 > \"", erlyRefSupportSams,"\" }; ",
      " if($0 ~ \"\tvv\"){print $0 > \"", erlyVarSupportSams,"\" } ",
    "}; r=$1}' OFS='\t' "
    #" }' OFS='\t' "
  )
  cmdString2 <- paste0(
    "vcf_read_support.py ",
    " -b ", lateBams,
    " -i ", vcf,
    if(quiet){" 2>/dev/null "},
    #" | awk '{ ",
    " | awk '{ if( !a[$1]++ ){",
      " if($0 ~ \"\tvr\"){print $0 > \"", lateRefSupportSams,"\" }; ",
      " if($0 ~ \"\tvv\"){print $0 > \"", lateVarSupportSams,"\" } ",
    "}; r=$1}' OFS='\t' "
    #" }' OFS='\t' "
  )

  dump<-mclapply(c(cmdString1,cmdString2),function(i){
    print(i)
    system(i)
  } , mc.cores=cores)

  erlyVarSupportBeds <- bed.sort(bam.2.bed(samtools.view(erlyVarSupportSams,genome.chrom.sizes=genomefile,includeHeader=T)))
  erlyRefSupportBeds <- bed.sort(bam.2.bed(samtools.view(erlyRefSupportSams,genome.chrom.sizes=genomefile,includeHeader=T)))
  lateVarSupportBeds <- bed.sort(bam.2.bed(samtools.view(lateVarSupportSams,genome.chrom.sizes=genomefile,includeHeader=T)))
  lateRefSupportBeds <- bed.sort(bam.2.bed(samtools.view(lateRefSupportSams,genome.chrom.sizes=genomefile,includeHeader=T)))

  varbed <- vcf.2.bed(vcf)

  # erlyVarSupportBgs <- unlist(mclapply(erlyVarSupportBeds,bed.coverage,windowfile=varbed,scalar=erlyFactors,mc.cores=cores))
  # erlyRefSupportBgs <- unlist(mclapply(erlyRefSupportBeds,bed.coverage,windowfile=varbed,scalar=erlyFactors,mc.cores=cores))
  # lateVarSupportBgs <- unlist(mclapply(lateVarSupportBeds,bed.coverage,windowfile=varbed,scalar=lateFactors,mc.cores=cores))
  # lateRefSupportBgs <- unlist(mclapply(lateRefSupportBeds,bed.coverage,windowfile=varbed,scalar=lateFactors,mc.cores=cores))

  cat("creating read density bedGraphs\n")
  erlyVarSupportBgs <- unlist(mclapply(1:numbams, function(x) bed.coverage ( erlyVarSupportBeds[x] , windowfile=varbed , scalar=erlyFactors[x] ),mc.cores=cores))
  erlyRefSupportBgs <- unlist(mclapply(1:numbams, function(x) bed.coverage ( erlyRefSupportBeds[x] , windowfile=varbed , scalar=erlyFactors[x] ),mc.cores=cores))
  lateVarSupportBgs <- unlist(mclapply(1:numbams, function(x) bed.coverage ( lateVarSupportBeds[x] , windowfile=varbed , scalar=lateFactors[x] ),mc.cores=cores))
  lateRefSupportBgs <- unlist(mclapply(1:numbams, function(x) bed.coverage ( lateRefSupportBeds[x] , windowfile=varbed , scalar=lateFactors[x] ),mc.cores=cores))

  cat("creating read density bigWigs\n")
  erlyVarSupportBws <- unlist(mclapply(erlyVarSupportBgs,bg.2.bw,genomefile=genomefile,mc.cores=cores))
  erlyRefSupportBws <- unlist(mclapply(erlyRefSupportBgs,bg.2.bw,genomefile=genomefile,mc.cores=cores))
  lateVarSupportBws <- unlist(mclapply(lateVarSupportBgs,bg.2.bw,genomefile=genomefile,mc.cores=cores))
  lateRefSupportBws <- unlist(mclapply(lateRefSupportBgs,bg.2.bw,genomefile=genomefile,mc.cores=cores))


  # define names of
  varDifBgsOutnames <- paste0(sampleNames,"_supportingVar_earlyMinusLate.bg")
  refDifBgsOutnames <- paste0(sampleNames,"_supportingRef_earlyMinusLate.bg")

  # define names of
  varLogBgsOutnames <- paste0(sampleNames,"_supportingVar_log2earlyOverLate.bg")
  refLogBgsOutnames <- paste0(sampleNames,"_supportingRef_log2earlyOverLate.bg")

  # define names of
  varMenBgsOutnames <- paste0(sampleNames,"_supportingVar_meanOfEandL.bg")
  refMenBgsOutnames <- paste0(sampleNames,"_supportingRef_meanOfEandL.bg")

  # make difference bedGraphs
  cat("creating read density difference bedGraphs\n")
  varDifBgs <- bg.ops(erlyVarSupportBgs,"difference",lateVarSupportBgs,outnames=varDifBgsOutnames)
  refDifBgs <- bg.ops(erlyRefSupportBgs,"difference",lateRefSupportBgs,outnames=refDifBgsOutnames)

  # make difference bedGraphs
  cat("creating mean read density bedGraphs\n")
  varMenBgs <- bg.ops(erlyVarSupportBgs,"mean",lateVarSupportBgs,outnames=varMenBgsOutnames)
  refMenBgs <- bg.ops(erlyRefSupportBgs,"mean",lateRefSupportBgs,outnames=refMenBgsOutnames)

  # make log2 bedGraphs
  cat("creating log2 ratio bedGraphs\n")
  varLogBgs <- bg.ops(erlyVarSupportBgs,"log2ratio",lateVarSupportBgs,outnames=varLogBgsOutnames)
  refLogBgs <- bg.ops(erlyRefSupportBgs,"log2ratio",lateRefSupportBgs,outnames=refLogBgsOutnames)

  # make big wigs
  cat("creating read density difference bigWigs\n")
  varDifBws <- unlist(lapply(varDifBgs,bg.2.bw,genomefile=genomefile))
  refDifBws <- unlist(lapply(refDifBgs,bg.2.bw,genomefile=genomefile))

  # make big wigs
  cat("creating log2 ratio bigWigs\n")
  varLogBws <- unlist(lapply(varLogBgs,bg.2.bw,genomefile=genomefile))
  refLogBws <- unlist(lapply(refLogBgs,bg.2.bw,genomefile=genomefile))

  # make big wigs
  cat("creating mean read density bigWigs\n")
  varMenBws <- unlist(lapply(varMenBgs,bg.2.bw,genomefile=genomefile))
  refMenBws <- unlist(lapply(refMenBgs,bg.2.bw,genomefile=genomefile))


  # erlyVarSupportScores <- read.bgs(erlyVarSupportBgs)
  # erlyRefSupportScores <- read.bgs(erlyRefSupportBgs)
  # lateVarSupportScores <- read.bgs(lateVarSupportBgs)
  # lateRefSupportScores <- read.bgs(lateRefSupportBgs)
  #
  # erlyVarSupportCounts<-unlist(lapply(erlyVarSupportSams,filelines))
  # erlyRefSupportCounts<-unlist(lapply(erlyRefSupportSams,filelines))
  # lateVarSupportCounts<-unlist(lapply(lateVarSupportSams,filelines))
  # lateRefSupportCounts<-unlist(lapply(lateRefSupportSams,filelines))
  #
  # erlyVarSupportRpm<-erlyVarSupportCounts*erlyFactors
  # erlyRefSupportRpm<-erlyRefSupportCounts*erlyFactors
  # lateVarSupportRpm<-lateVarSupportCounts*lateFactors
  # lateRefSupportRpm<-lateRefSupportCounts*lateFactors
  #
  # VarsupportLog2<-log2(erlyVarSupportRpm/lateVarSupportRpm)
  # RefSupportLog2<-log2(erlyRefSupportRpm/lateRefSupportRpm)

  # out <- data.frame(
  #   sampleNames,
  #   erlyBams,
  #   lateBams,
  #   erlyCounts,
  #   lateCounts,
  #   erlyVarSupportCounts,
  #   lateVarSupportCounts,
  #   erlyRefSupportCounts,
  #   lateRefSupportCounts,
  #   erlyVarSupportRpm,
  #   lateVarSupportRpm,
  #   erlyRefSupportRpm,
  #   lateRefSupportRpm,
  #   VarsupportLog2,
  #   RefSupportLog2
  # )
  #
  #
  # #out<-(list(SupportSams,RefSupportSams))
  # #names(out) <- c("ref","var")
  # return(out)
}
