#!/bin/bash

docker run --rm -it -v `pwd`:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib  -w /data trinityctat/ctat_lr_fusion:latest $*

