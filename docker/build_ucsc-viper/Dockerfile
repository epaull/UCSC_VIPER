FROM epaull/r-latest

# Building Scipy may need more than 1G RAM to complete succssfully.
# Some minimal cloud VPS servers are not sufficient.

RUN apt-get -y update && apt-get install -y git
RUN apt-get -y update && apt-get install -y make


# update the sources list 
#RUN apt-get -y update && apt-get install -y r-base r-base-dev

RUN apt-get -y install libcurl4-openssl-dev
WORKDIR /home
#setup R configs
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('mygene')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('Biobase'); install.packages('Biobase')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('mixtools'); install.packages('mixtools')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('viper'); install.packages('viper')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('GenomicFeatures')"
RUN Rscript -e "install.packages('getopt')"

RUN rm -rf UCSC_VIPER
RUN git clone http://github.com/epaull/UCSC_VIPER.git

WORKDIR /data
ENTRYPOINT ["sh", "/home/UCSC_VIPER/docker/build_ucsc-viper/viper_wrapper.sh"]
