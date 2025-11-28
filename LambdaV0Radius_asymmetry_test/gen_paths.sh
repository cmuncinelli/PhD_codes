#!/bin/bash

RUN_NUMBER=$1
MAX_INDEX=$2

BASE="/alice/cern.ch/user/a/alihyperloop/jobs/0309"

for i in $(seq -f "%03g" 1 $MAX_INDEX); do
    echo "${BASE}/hy_${RUN_NUMBER}/AOD/${i}/AO2D.root"
done

# Runs like this:
# ./gen_paths.sh 3090064 39