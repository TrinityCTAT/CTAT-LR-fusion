#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

cachebustdate=$(date +%s)

docker build -t trinityctat/ctat_lr_fusion:${VERSION} --build-arg CACHEBUST=${cachebustdate} .
docker build -t trinityctat/ctat_lr_fusion:latest --build-arg CACHEBUST=${cachebustdate} .

