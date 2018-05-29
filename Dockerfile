FROM ubuntu:16.04

MAINTAINER Kyle Hernandez <kmhernan@uchicago.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y --force-yes \
        build-essential \
        python3.5 \
        python3.5-dev \
        python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt

## Install python package
RUN mkdir /opt/gdc-rnaseq-tool
WORKDIR /opt/gdc-rnaseq-tool
ADD utils /opt/gdc-rnaseq-tool/
ADD LICENSE /opt/gdc-rnaseq-tool/

WORKDIR /opt
