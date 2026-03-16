#!/usr/bin/env bash

###############################################################################
# Script name:
#   convert_list_of_paths.sh
#
# Purpose:
#   Scans all XML files inside the aod_collections/ subfolder of a given
#   working directory, extracts all Logical File Name (lfn) paths found in
#   <file ... /> entries, and writes them one per line into:
#
#     <WORK_DIR>/DownloadListFromAlien.txt
#
#   If the output file already exists it is overwritten.
#
# Arguments:
#   $1 : WORK_DIR -- absolute path to the per-wagon working directory.
#        This directory must contain an aod_collections/ subfolder with
#        one or more XML files from MonALISA.
#
# Example:
#   ./convert_list_of_paths.sh /home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC_min3ITS
#
# Output example:
#   /alice/data/2025/LHC25ae/564445/apass2/1520/.../AO2D.root
#
###############################################################################

# Stop the script if any command fails
set -e

###############################################################################
# Read and validate the WORK_DIR argument
###############################################################################
WORK_DIR="$1"

if [ -z "$WORK_DIR" ]; then
  echo "Usage: $0 <WORK_DIR>"
  echo "  WORK_DIR : path to the per-wagon working directory"
  exit 1
fi

if [ ! -d "$WORK_DIR" ]; then
  echo "Error: WORK_DIR does not exist: ${WORK_DIR}"
  exit 1
fi

###############################################################################
# Derive paths from WORK_DIR
###############################################################################
INPUT_DIR="${WORK_DIR}/aod_collections"
OUTPUT_FILE="${WORK_DIR}/DownloadListFromAlien.txt"

if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: aod_collections/ not found inside: ${WORK_DIR}"
  exit 1
fi

###############################################################################
# Clear any previous output file to avoid accidental appending
###############################################################################
rm -f "${OUTPUT_FILE}"

###############################################################################
# Glob for XML files; nullglob prevents a literal "*.xml" match on empty dirs
###############################################################################
shopt -s nullglob
XML_FILES=("${INPUT_DIR}"/*.xml)

if [ ${#XML_FILES[@]} -eq 0 ]; then
  echo "No XML files found inside ${INPUT_DIR}."
  exit 0
fi

###############################################################################
# Process each XML file
###############################################################################
for XML_FILE in "${XML_FILES[@]}"; do

  echo "Processing ${XML_FILE} ..."

  # Extract the value of the lfn attribute from lines containing lfn="
  #
  # Pipeline:
  #   grep 'lfn="'                   -- select only lines with the lfn attribute
  #   sed 's/.*lfn="\([^"]*\)".*/\1/'  -- capture the path between the quotes
  #
  # Result: one clean AliEN path per line, e.g.
  #   /alice/data/2025/LHC25ae/564445/apass2/1520/.../AO2D.root
  #
  grep 'lfn="' "${XML_FILE}" \
    | sed 's/.*lfn="\([^"]*\)".*/\1/' \
    >> "${OUTPUT_FILE}"

done

echo "Done. LFN list written to ${OUTPUT_FILE}"
