FROM ubuntu:14.04

# Building Scipy may need more than 1G RAM to complete succssfully.
# Some minimal cloud VPS servers are not sufficient.

RUN apt-get -y update && apt-get install -y git
RUN apt-get -y update && apt-get install -y make
RUN apt-get -y update && apt-get install -y python python-pip python-dev python-numpy 
RUN apt-get -y update && apt-get install -y python-setuptools python-dateutil 
RUN apt-get -y update && apt-get install -y python-scipy

## Install the r-base stuff
RUN pip install pandas
RUN apt-get install -y r-recommended
RUN apt-get -y update && apt-get install -y python-rpy2


#RUN git clone https://github.com/ucscCancer/pathmark-scripts.git
#RUN cd pathmark-scripts && make
#RUN cd pathmark-scripts && git clone https://github.com/ucscCancer/paradigm-scripts.git
#RUN cd pathmark-scripts/paradigm-scripts && make
#RUN cd pathmark-scripts/paradigm-scripts && source init.csh
#RUN cd pathmark-scripts && source init.csh
