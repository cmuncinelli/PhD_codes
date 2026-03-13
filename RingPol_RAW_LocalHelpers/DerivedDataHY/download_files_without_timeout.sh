#!/usr/bin/env bash

###############################################################################
# Script name:
#   download_files_without_timeout.sh
#
# Purpose:
#   Generic AliEN file downloader. Downloads every LFN listed in a given
#   input file into a given output directory, naming files sequentially as:
#
#     <FILE_PREFIX>_1.root
#     <FILE_PREFIX>_2.root
#     ...
#
#   Files that already exist in the output directory are skipped, making
#   interrupted downloads safely resumable.
#
#   Timeout is intentionally disabled: Hyperloop-derived AODs and
#   AnalysisResults files can be several GB each, and alien_cp manages
#   its own transfer retries internally.
#
# Arguments:
#   $1 : LIST_FILE   -- full path to a plain-text file with one AliEN LFN
#                       per line (as produced by convert_alien_paths.sh)
#   $2 : OUTPUT_DIR  -- full path to the directory where files will be saved
#   $3 : FILE_PREFIX -- prefix for the output filenames (e.g. "AO2D" or
#                       "AnalysisResults")
#
# Example (AO2D files):
#   ./download_files_without_timeout.sh \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2_small/ITSandTPC/DownloadListFromAlien.txt \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2_small/ITSandTPC/AO2Ds \
#       AO2D
#
# Example (AnalysisResults files):
#   ./download_files_without_timeout.sh \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2_small/ITSandTPC/DownloadListAnalysisResults.txt \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2_small/ITSandTPC/temp_analysis_results \
#       AnalysisResults
#
###############################################################################

###############################################################################
# Helper: strip any trailing slashes from a path
###############################################################################
strip_trailing_slashes() {
  local path="$1"
  path="${path%/}"
  echo "$path"
}

###############################################################################
# Read and validate arguments
###############################################################################
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  echo "Usage: $0 <LIST_FILE> <OUTPUT_DIR> <FILE_PREFIX>"
  echo ""
  echo "  LIST_FILE   : full path to the AliEN LFN list"
  echo "  OUTPUT_DIR  : full path to the destination directory"
  echo "  FILE_PREFIX : filename prefix, e.g. AO2D or AnalysisResults"
  exit 1
fi

LIST_FILE="$(strip_trailing_slashes "$1")"
OUTPUT_DIR="$(strip_trailing_slashes "$2")"
FILE_PREFIX="$3"

if [ ! -f "$LIST_FILE" ]; then
  echo "Error: LIST_FILE not found: ${LIST_FILE}"
  exit 1
fi

# Create output directory if it does not yet exist
mkdir -p "$OUTPUT_DIR"

###############################################################################
# Main download loop
###############################################################################
COUNTER=1
MAX=$(wc -l < "$LIST_FILE")

while IFS= read -r FILE || [ -n "$FILE" ]; do

  echo "[FILE #${COUNTER}/${MAX}: ${FILE}]"

  OUTFILE="${OUTPUT_DIR}/${FILE_PREFIX}_${COUNTER}.root"

  if [ -f "$OUTFILE" ]; then
    echo "  ${FILE_PREFIX}_${COUNTER}.root already exists -- skipping."
  else
    echo "  Copying: alien:${FILE}  ->  ${OUTFILE}"

    # alien_cp options:
    #   -t 6 : use up to 6 parallel transfer threads per file
    #
    # To re-enable a hard timeout (e.g. for smaller files), replace with:
    #   timeout 120s alien_cp -t 6 "alien:${FILE}" "file:${OUTFILE}"
    # and handle RC=124 accordingly.
    alien_cp -t 6 "alien:${FILE}" "file:${OUTFILE}"
  fi

  (( ++COUNTER ))

done < "$LIST_FILE"

echo ""
echo "Download complete. Files are in: ${OUTPUT_DIR}"
