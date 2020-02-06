#!/bin/bash
# Builds the Dockerfile and tags with current git commit. Only
# relevant for internal GDC building.

docker_build () {
    command docker build \
        --build-arg http_proxy=$http_proxy \
        --build-arg https_proxy=$https_proxy \
        --build-arg no_proxy=$no_proxy \
        "$@"
}

# tag
quay="quay.io/ncigdc/gdc-rnaseq-tool"
version=$(git log --first-parent --max-count=1 --format=format:%H)
imagetag="${quay}:${version}"

echo "Building tag: $imagetag"
docker_build -t $imagetag .
