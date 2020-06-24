#!/bin/sh
../nextflow run ./main.nf \
-profile gls \
-c /home/combizkhozoie/nf-core-scflow/conf/gls.config \
-params-file /home/combizkhozoie/scflowexamplegcp/conf/nfx-params.json \
-with-tower \
-resume
