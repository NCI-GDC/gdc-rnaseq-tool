ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.8-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /gdc_rnaseq_tools

WORKDIR /gdc_rnaseq_tools

RUN pip install tox && tox -e build

FROM ${REGISTRY}/python3.8:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="gdc_rnaseq_tools" \
      org.opencontainers.image.description="Utility scripts for GDC RNA-seq workflows. The docker file also installs Trimmomatic and fqvendorfail." \
      org.opencontainers.image.source="https://github.com/NCI-GDC/gdc-rnaseq-tool" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /gdc_rnaseq_tools/dist/*.whl /gdc_rnaseq_tools/
COPY requirements.txt /gdc_rnaseq_tools/

WORKDIR /gdc_rnaseq_tools

RUN pip install --no-deps -r requirements.txt \
	&& pip install --no-deps *.whl \
	&& rm -f *.whl requirements.txt

CMD ["gdc_rnaseq_tools --help"]
