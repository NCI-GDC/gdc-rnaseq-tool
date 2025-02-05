ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.6-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /gdc_rnaseq_tool

WORKDIR /gdc_rnaseq_tool

RUN pip install tox && tox -e build

FROM ${REGISTRY}/python3.6:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="gdc_rnaseq_tool" \
      org.opencontainers.image.description="Utility scripts for GDC RNA-seq workflows. The docker file also installs Trimmomatic and fqvendorfail." \
      org.opencontainers.image.source="https://github.com/NCI-GDC/gdc-rnaseq-tool" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /gdc_rnaseq_tool/dist/*.whl /gdc_rnaseq_tool/
COPY requirements.txt /gdc_rnaseq_tool/

WORKDIR /gdc_rnaseq_tool

RUN pip install --no-deps -r requirements.txt \
	&& pip install --no-deps *.whl \
	&& rm -f *.whl requirements.txt

USER app

CMD ["--help"]
