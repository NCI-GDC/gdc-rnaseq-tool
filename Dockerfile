FROM quay.io/ncigdc/bio-python:3.6 as builder

COPY ./ /opt

WORKDIR /opt

RUN python -m pip install tox && tox

# tox step builds sdist

FROM quay.io/ncigdc/bio-python:3.6

COPY --from=builder /opt/dist/*.tar.gz /opt
COPY ./requirements.txt /opt/requirements.txt

WORKDIR /opt

# Install package from sdist
RUN pip install -r requirements.txt \
	&& pip install *.tar.gz \
	&& rm -rf *.tar.gz requirements.txt

ENTRYPOINT ["gdc_rnaseq_tools"]

CMD ["--help"]
