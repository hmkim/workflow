FROM python:2.7
MAINTAINER john.didion@nih.gov

RUN apt-get update && apt-get -y install libcurl4-nss-dev

ADD https://github.com/arq5x/bedtools2/releases/download/v2.23.0/bedtools-2.25.0.tar.gz bedtools-2.25.0.tar.gz
RUN tar zxf bedtools-2.25.0.tar.gz
RUN cd bedtools2 && make && make install
RUN rm -rf bedtools*