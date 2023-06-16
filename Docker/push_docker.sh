#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/minimap2fusion:${VERSION}
docker push trinityctat/minimap2fusion:latest

