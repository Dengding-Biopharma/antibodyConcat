FROM continuumio/miniconda3

ARG http_proxy
ARG https_proxy

WORKDIR /home

RUN http_proxy=$http_proxy https_proxy=$https_proxy conda install -c bioconda rapsearch

RUN apt-get update
RUN apt-get install -y unzip g++

# install dependencies
RUN mkdir -p /install/

WORKDIR /install
RUN http_proxy=$http_proxy https_proxy=$https_proxy wget https://github.com/lkytal/PredFull/archive/refs/heads/master.zip -O PredFull.zip
RUN unzip PredFull.zip
WORKDIR /install/PredFull-master
RUN pip3 install pandas==1.5.2 tensorflow==2.11.0 pyteomics==4.5.6 lxml==4.9.1 -i https://pypi.doubanio.com/simple/

WORKDIR /install
RUN http_proxy=$http_proxy https_proxy=$https_proxy wget https://github.com/COL-IU/msSLASH/archive/refs/heads/master.zip -O msSLASH.zip
RUN unzip msSLASH.zip
WORKDIR /install/msSLASH-master
RUN bash install.sh

WORKDIR /home
ADD *.py .

# COPY requirements .
# RUN pip install -r requirements.txt -i https://pypi.doubanio.com/simple/
