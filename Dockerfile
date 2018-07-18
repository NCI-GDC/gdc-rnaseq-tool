FROM ubuntu:16.04

MAINTAINER Kyle Hernandez <kmhernan@uchicago.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y --force-yes \
        build-essential \
        python3.5 \
        python3.5-dev \
        python3-pip \
        wget \
        unzip \
        openjdk-8-jre-headless \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt

## Install trimmomatic to /opt/Trimmomatic-0.38/trimmomatic-0.38.jar
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip \
    && unzip Trimmomatic-0.38.zip \
    && rm Trimmomatic-0.38.zip

## Install python package
WORKDIR /opt
RUN mkdir /opt/gdc-rnaseq-tool
WORKDIR /opt/gdc-rnaseq-tool
ADD utils /opt/gdc-rnaseq-tool/
ADD LICENSE /opt/gdc-rnaseq-tool/

## Install fqvendorfail
WORKDIR /opt
RUN git clone git@github.com:kmhernan/fqvendorfail.git \
    && cd fqvendorfail \
    && make
