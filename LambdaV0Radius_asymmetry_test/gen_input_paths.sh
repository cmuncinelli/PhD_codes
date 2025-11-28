#!/bin/bash

MAX_INDEX=$1
BASE="file:/storage2/cicero/asymmetric_rapidity_V0tests_LHC23_pass4_Thin_small"
# Added "file:" as a prefix because O2 was interpreting the input_data.txt list as a list of LFN entries at the GRID!

for i in $(seq 1 $MAX_INDEX); do
    echo "${BASE}/AO2D_${i}.root"
done

# Runs like ./gen_input_paths.sh 30