FROM ubuntu:16.04

MAINTAINER Kyle Hernandez <kmhernan@uchicago.edu>

RUN DEBIAN_FRONTEND=noninteractive apt-get -y update && apt-get install -y --force-yes \
        build-essential \
        python3.5 \
        python3.5-dev \
        python3-pip \
        wget \
        unzip \
        openjdk-8-jre-headless \
        git \
        zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install trimmomatic to /opt/Trimmomatic-0.38/trimmomatic-0.38.jar
RUN cd /opt \
    && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip \
    && unzip Trimmomatic-0.38.zip \
    && rm Trimmomatic-0.38.zip

## Install python package
RUN mkdir /opt/gdc-rnaseq-tool \
    && cd /opt/gdc-rnaseq-tool
ADD utils /opt/gdc-rnaseq-tool/
ADD LICENSE /opt/gdc-rnaseq-tool/

## Install fqvendorfail
RUN cd /opt \
    && git clone https://github.com/kmhernan/fqvendorfail.git \
    && cd fqvendorfail \
    && git checkout 29ca20856a33393cd57a022dfd5687cea18332f7 \
    && make

WORKDIR /opt
