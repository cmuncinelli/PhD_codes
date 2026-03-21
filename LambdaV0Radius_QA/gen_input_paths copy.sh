#!/bin/bash

MAX_INDEX=$1
BASE="file:/storage2/cicero/asymmetric_rapidity_V0tests_LF_LHC23k4f_pass4_Strangeness"
# Added "file:" as a prefix because O2 was interpreting the input_data.txt list as a list of LFN entries at the GRID!

OUTPUT="input_data_storage.txt"

# Clear output file before writing
> "$OUTPUT"

for i in $(seq 1 "$MAX_INDEX"); do
    echo "${BASE}/AO2D_${i}.root" >> "$OUTPUT"
done

# Runs like ./gen_input_paths.sh 30