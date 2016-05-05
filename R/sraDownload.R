sraDownload <- function(srp){
  system("wget -O- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP057023' | tr ',' '\t'") 
}
