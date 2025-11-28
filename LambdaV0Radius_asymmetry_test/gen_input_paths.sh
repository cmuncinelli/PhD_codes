#!/bin/bash

MAX_INDEX=$1
BASE="/storage2/cicero/asymmetric_rapidity_V0tests_LHC23_pass4_Thin_small"

for i in $(seq 1 $MAX_INDEX); do
    echo "${BASE}/AO2D_${i}.root"
done

# Runs like ./gen_input_paths.sh 30