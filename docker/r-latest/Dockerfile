FROM ubuntu:14.04

# add newer sources
RUN sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
# CRAN archives key
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
# add key to APT
RUN gpg -a --export E084DAB9 | sudo apt-key add -
#RUN sudo apt-get update

# install latest R
RUN apt-get -y update && apt-get install -y r-base r-base-dev

