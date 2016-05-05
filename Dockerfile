FROM rocker/hadleyverse
MAINTAINER "Daniel Vera" vera@genomics.fsu.edu

ENV SWDIR=/opt
ENV PATH /opt/bowtie:/opt/samtools-1.2:/opt/bedtools/bin:/opt/bcftools-1.2:/opt/kent:$PATH

RUN apt-get update
#RUN apt-get upgrade -y
RUN apt-get install -y python-dev python-pip man openssh-client
RUN mkdir -p $SWDIR
RUN git clone http://github.com/arq5x/bedtools2.git $SWDIR/bedtools &&\
cd $SWDIR/bedtools &&\
make
#RUN mkdir $SWDIR/kent &&\
#cd $SWDIR/kent &&\
#wget -nv -r -np -e robots=off -R 'html' http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/ &&\
#mv hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/* . &&\
#rm index* &&\
#rm FOOTER &&\
#rm -rf hgdownload.cse.ucsc.edu &&\
#chmod +x *
RUN git clone https://github.com/BenLangmead/bowtie.git $SWDIR/bowtie &&\
cd $SWDIR/bowtie &&\
make
RUN cd $SWDIR &&\
wget -O- https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 | bunzip2 -c | tar -xv &&\
cd samtools-1.2 &&\
make
RUN cd $SWDIR &&\
wget -O- https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2 | bunzip2 -c | tar -xv &&\
cd bcftools-1.2 &&\
make
RUN Rscript -e "devtools::install_github('dvera/rage')"
#RUN pip install numpy
#RUN pip install macs2

#echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
#Rscript -e "install.packages('devtools')
#Rscript -e "install.packages('gtools')
#Rscript -e "install.packages('gplots')



# bgzip
# tabix
# bowtie2
# tophat
# hisat
# cufflinks
# bwa
# seqtk
# edgeR
# DESeq2
# fastqc
# bfc
# fastx_toolkit
# gunzip?
# sshfs
