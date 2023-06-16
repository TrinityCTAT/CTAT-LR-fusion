#!/bin/bash

docker run --rm -it -v `pwd`:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib  trinityctat/gmapfusion:latest $*

