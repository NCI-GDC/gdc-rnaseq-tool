FROM quay.io/ncigdc/bio-python:3.6

ENV BINARY=gdc_rnaseq_tools

RUN apt-get update \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*

COPY ./dist/ /opt

WORKDIR /opt

RUN make init-pip \
  && ln -s /opt/bin/${BINARY} /bin/${BINARY} \
  && chmod +x /bin/${BINARY}

ENTRYPOINT ["/bin/gdc_rnaseq_tools"]

CMD ["--help"]
