#!/usr/bin/env bash

###############################################################################
# Script name:
#   convert_list_of_paths.sh
#
# Purpose:
#   This script scans one or more XML files inside the ./aod_collections/ directory
#   and extracts all Logical File Name (lfn) paths found in <file ... /> entries.
#   These paths correspond to AO2D.root files stored in the ALICE grid.
#
#   The script collects all extracted lfn paths and writes them, one per line,
#   into a plain text file named:
#
#     DownloadListFromAlien.txt
#
#   The output file is created in the same directory where the script is run.
#   If the output file already exists, it is overwritten.
#
# Expected input:
#   - One or more XML files in the current directory.
#   - Typical examples are:
#       aod_collection.xml
#       aod_collection_1520.xml
#
# Output:
#   - A text file named DownloadListFromAlien.txt containing full ALICE paths,
#     for example:
#       /alice/data/2025/LHC25ae/564445/apass2/1520/.../AO2D.root
#
# How to use: ./convert_list_of_paths.sh
#
###############################################################################

# Stop the script if any command fails
set -e

# Directory containing XML collections
INPUT_DIR="aod_collections"

# Check if input directory exists
if [ ! -d "${INPUT_DIR}" ]; then
  echo "Error: Directory '${INPUT_DIR}' does not exist."
  exit 1
fi

# Name of the output file that will contain the list of LFNs
OUTPUT_FILE="DownloadListFromAlien.txt"

# Remove the output file if it already exists
# This avoids appending to an old file by accident
rm -f "${OUTPUT_FILE}"

# Enable nullglob so that *.xml expands to empty if no match
shopt -s nullglob

XML_FILES=("${INPUT_DIR}"/*.xml)

# Check if no XML files found
if [ ${#XML_FILES[@]} -eq 0 ]; then
  echo "No XML files found inside ${INPUT_DIR}."
  exit 0
fi

# Loop over XML files
for XML_FILE in "${XML_FILES[@]}"; do

  # If there are no XML files, the glob will expand to "*.xml"
  # This check avoids trying to process a non-existing file
  if [ ! -f "${XML_FILE}" ]; then
    continue
  fi

  # Inform the user which file is being processed
  echo "Processing ${XML_FILE} ..."

  # Extract the value of the lfn attribute from lines that contain lfn="
  #
  # Explanation of the pipeline:
  #
  # 1. grep 'lfn="' :
  #    Selects only lines that contain the lfn attribute
  #
  # 2. sed 's/.*lfn="\([^"]*\)".*/\1/' :
  #    Uses a regular expression to:
  #    - Remove everything before lfn="
  #    - Capture everything up to the next double quote
  #    - Remove everything after that
  #
  # The final result is a clean path like:
  # /alice/data/2025/LHC25ae/564445/apass2/1520/.../AO2D.root
  #
  grep 'lfn="' "${XML_FILE}" \
    | sed 's/.*lfn="\([^"]*\)".*/\1/' \
    >> "${OUTPUT_FILE}"

done

# Final message to confirm successful completion
echo "Done."
echo "LFN list written to ${OUTPUT_FILE}"
