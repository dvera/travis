bed.avglengths <-
function( databedfile , targetfile , windowsize=1000000 , genome=TRUE , fraglenoutfile="default" , numfragoutfile="default"){

	library(tools)

	# make sure input is only 1 file
	if(length(databedfile) > 1){stop("bed.parselengths can only take 1 file")}

	# prepare output file information
   ext<-file_ext(databedfile)
	fragname<-basename(removeext(databedfile))
	if(fraglenoutfile=="default"){
      fraglenoutfile = paste(fragname,"_",windowsize,"bp_win_frag_lens.bed",sep='')
   }
	if(numfragoutfile=="default"){
      numfragoutfile = paste(fragname,"_",windowsize,"bp_win_frag_counts.bed",sep='')
   }

   # if a target bed file is included, use it to isolate regions
   if(genome){
      # make windows from the chrom.sizes file
      cat("Calculating windows for genome\n")
      bed.makewindows(targetfile,windowsize=windowsize,genome=TRUE)
   }
   else{
      # make windows from the seqcap targets file
      cat("Calculating windows for target regions\n")
      bed.makewindows(targetfile,windowsize=windowsize)
   }
   
   # figure out name of windowed target file
   windowedtargets<-paste(basename(removeext(targetfile)),"_win",windowsize,".bed",sep="")
   
   # calculate average fragment length and number of fragments
   cat("Calculating fragment statistics\n")
   system(paste("awk '{print $1,$2,$3,$3-$2}' OFS='\t' ",databedfile,">",paste(fraglenoutfile,".all",sep='')))

   system(paste("bedtools map -a",windowedtargets,"-b",paste(fraglenoutfile,".all",sep=''),"-c 4 -null \"NA\" -o mean | grep -v NA >",fraglenoutfile))
   system(paste("bedtools map -a",windowedtargets,"-b",paste(fraglenoutfile,".all",sep=''),"-c 4 -null \"NA\" -o count | grep -v NA >",numfragoutfile))
   
   # clean up temp file with all fragment lengths (MAYBE)
   #system(paste("rm",paste(fraglenoutfile,".all",sep='')))
   
}
