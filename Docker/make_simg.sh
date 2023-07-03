#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build ctat_lr_fusion.v${VERSION}.simg docker://trinityctat/ctat_lr_fusion:$VERSION

singularity exec -e ctat_lr_fusion.v${VERSION}.simg ctat-LR-fusion

ln -sf  ctat_lr_fusion.v${VERSION}.simg  ctat_lr_fusion.simg

