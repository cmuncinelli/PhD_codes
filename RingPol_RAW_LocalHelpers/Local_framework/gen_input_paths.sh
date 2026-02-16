#!/usr/bin/env bash

###############################################################################
# Script name:
#   gen_input_paths.sh
#
# Purpose:
#   This script generates a local input file list for O2 workflows.
#   It assumes that AO2D files have already been downloaded to a local
#   storage directory and are named sequentially as:
#
#     AO2D_1.root
#     AO2D_2.root
#     ...
#     AO2D_N.root
#
#   The script builds full file paths using a user-provided storage path
#   and writes them into a text file named:
#
#     input_data_storage.txt
#
# Arguments:
#   $1 : Base storage path where AO2D files are located
#   $2 : Maximum index (number of AO2D files)
#
# Example usage:
#   ./gen_input_paths.sh /storage2/cicero/my_analysis 134
#
###############################################################################

# Stop the script if any command fails
set -e

# Read arguments
STORAGE_PATH="$1"
MAX_INDEX="$2"

# Safety check: make sure both arguments are provided
if [ -z "$STORAGE_PATH" ] || [ -z "$MAX_INDEX" ]; then
  echo "Usage: $0 <storage_path> <max_index>"
  exit 1
fi

# Build BASE path
# Added "file:" prefix because O2 interprets plain paths as GRID LFNs
BASE="file:${STORAGE_PATH}"

# Name of the output file
OUTPUT="input_data_storage.txt"

# Clear output file before writing
> "$OUTPUT"

# Generate the list of input paths
for i in $(seq 1 "$MAX_INDEX"); do
  echo "${BASE}/AO2D_${i}.root" >> "$OUTPUT"
done

# Inform the user
echo "Generated ${MAX_INDEX} input paths in ${OUTPUT}"