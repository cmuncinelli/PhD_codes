#!/usr/bin/env bash

###############################################################################
# Script name:
#   gen_input_paths.sh
#
# Purpose:
#   Generates a local input file list for O2 workflows.
#   Assumes AO2D files have already been (or will be) downloaded to:
#
#     <WORK_DIR>/AO2Ds/
#
#   and are named sequentially:
#     AO2D_1.root, AO2D_2.root, ..., AO2D_N.root
#
#   The output is written to:
#     <WORK_DIR>/input_data_storage.txt
#
#   Each line is prefixed with "file:" so O2 interprets them as local paths
#   rather than GRID LFNs.
#
# Arguments:
#   $1 : WORK_DIR  -- absolute path to the per-wagon working directory
#   $2 : MAX_INDEX -- number of AO2D files (i.e. lines in DownloadListFromAlien.txt)
#
# Example:
#   ./gen_input_paths.sh /home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC_min3ITS 134
#
###############################################################################

# Stop the script if any command fails
set -e

###############################################################################
# Read and validate arguments
###############################################################################
WORK_DIR="$1"
MAX_INDEX="$2"

if [ -z "$WORK_DIR" ] || [ -z "$MAX_INDEX" ]; then
  echo "Usage: $0 <WORK_DIR> <MAX_INDEX>"
  echo "  WORK_DIR  : path to the per-wagon working directory"
  echo "  MAX_INDEX : number of AO2D files to list"
  exit 1
fi

if [ ! -d "$WORK_DIR" ]; then
  echo "Error: WORK_DIR does not exist: ${WORK_DIR}"
  exit 1
fi

###############################################################################
# Derive paths
#
# AO2D files live in the AO2Ds/ subfolder.
# The "file:" prefix tells O2 these are local filesystem paths, not grid LFNs.
###############################################################################
AOD_DIR="${WORK_DIR}/AO2Ds"
OUTPUT="${WORK_DIR}/input_data_storage.txt"

###############################################################################
# Write the file list
###############################################################################
> "$OUTPUT"

for i in $(seq 1 "$MAX_INDEX"); do
  echo "file:${AOD_DIR}/AO2D_${i}.root" >> "$OUTPUT"
done

echo "Generated ${MAX_INDEX} input paths in ${OUTPUT}"
