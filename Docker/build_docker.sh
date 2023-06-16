#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

cachebustdate=$(date +%s)

docker build -t trinityctat/minimap2fusion:${VERSION} --build-arg CACHEBUST=${cachebustdate} .
docker build -t trinityctat/minimap2fusion:latest --build-arg CACHEBUST=${cachebustdate} .

