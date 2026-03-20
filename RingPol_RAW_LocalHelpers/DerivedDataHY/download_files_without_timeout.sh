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
# Files that already exist in the output directory are skipped, making
#   interrupted downloads safely resumable.
#
#   A machine-readable stats file is optionally written to $4. It contains
#   KEY=VALUE lines consumed by download_hyperloop.sh to build the summary log.
#
# Arguments:
#   $1 : LIST_FILE   -- full path to a plain-text file with one AliEN LFN
#                       per line (as produced by convert_alien_paths.sh)
#   $2 : OUTPUT_DIR  -- full path to the destination directory
#   $3 : FILE_PREFIX -- prefix for the output filenames (e.g. "AO2D" or
#                       "AnalysisResults")
#   $4 : STATS_FILE  -- (optional) path where KEY=VALUE download stats are
#                       written for consumption by the orchestrator
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
  echo "Usage: $0 <LIST_FILE> <OUTPUT_DIR> <FILE_PREFIX> [STATS_FILE]"
  echo ""
  echo "  LIST_FILE   : full path to the AliEN LFN list"
  echo "  OUTPUT_DIR  : full path to the destination directory"
  echo "  FILE_PREFIX : filename prefix, e.g. AO2D or AnalysisResults"
  echo "  STATS_FILE  : (optional) path for machine-readable KEY=VALUE stats"
  exit 1
fi

LIST_FILE="$(strip_trailing_slashes "$1")"
OUTPUT_DIR="$(strip_trailing_slashes "$2")"
FILE_PREFIX="$3"
STATS_FILE="${4:-}"

if [ ! -f "$LIST_FILE" ]; then
  echo "Error: LIST_FILE not found: ${LIST_FILE}"
  exit 1
fi

# Create output directory if it does not yet exist
mkdir -p "$OUTPUT_DIR"

###############################################################################
# Counters
###############################################################################
TOTAL=0
SKIPPED=0
DOWNLOADED=0
FAILED=0

###############################################################################
# Main download loop
###############################################################################
MAX=$(wc -l < "$LIST_FILE")

while IFS= read -r FILE || [ -n "$FILE" ]; do

  [ -z "$FILE" ] && continue
  (( TOTAL += 1 ))

  echo "[${FILE_PREFIX} #${TOTAL}/${MAX}: ${FILE}]"

  OUTFILE="${OUTPUT_DIR}/${FILE_PREFIX}_${TOTAL}.root"

  # Soft safety net: skip any LFN that contains a known non-directory
  # entry name (e.g. aod_collection.xml). These should have been filtered
  # by convert_alien_paths.sh, but this guard prevents a bad entry in an
  # existing DownloadList file from causing a failed alien_cp.
  # If we don't do this soft skip, I will need to redownload all derived
  # data files (a real pain!), so kept the old lists to avoid index mis-
  # matching. We are now avoiding errors like this, whilst still preserving
  # the names of all files that would come after the attempted "AO2D_66.root":
  #   [AO2D #66/105: /alice/cern.ch/user/a/alihyperloop/jobs/0441/hy_4410967/AOD/aod_collection.xml/AO2D.root]
  #   Copying: alien:/alice/cern.ch/user/a/alihyperloop/jobs/0441/hy_4410967/AOD/aod_collection.xml/AO2D.root  ->  /home/users/cicerodm/RingPol/LHC25ae_pass2_small/ITSandTPC_min2ITS/AO2Ds/AO2D_66.root
  #   Specified source /alice/cern.ch/user/a/alihyperloop/jobs/0441/hy_4410967/AOD/aod_collection.xml/AO2D.root not found!
  #   [WARN] alien_cp failed (exit 2).

  if echo "$FILE" | grep -q "aod_collection.xml"; then
    echo "  Skipping non-directory entry: ${FILE}"
    (( SKIPPED += 1 ))
    continue
  fi

  if [ -f "$OUTFILE" ]; then
    echo "  Already exists -- skipping."
    (( SKIPPED += 1 ))
    continue
  fi

  echo "  Copying: alien:${FILE}  ->  ${OUTFILE}"

  # alien_cp options:
  #   -t 6 : use up to 6 parallel transfer threads per file
  #
  # Exit code is captured explicitly; set -e is NOT used here so that a
  # single failed file does not abort the entire download run.
  alien_cp -t 6 "alien:${FILE}" "file:${OUTFILE}"
  CP_EXIT=$?

  if [ $CP_EXIT -ne 0 ]; then
    echo "  [WARN] alien_cp failed (exit ${CP_EXIT})."
    (( FAILED += 1 ))
  else
    (( DOWNLOADED += 1 ))
  fi

done < "$LIST_FILE"

###############################################################################
# Write stats file if requested
###############################################################################
if [ -n "$STATS_FILE" ]; then
  {
    echo "TOTAL=${TOTAL}"
    echo "SKIPPED=${SKIPPED}"
    echo "DOWNLOADED=${DOWNLOADED}"
    echo "FAILED=${FAILED}"
  } > "$STATS_FILE"
fi

echo ""
echo "Download complete. Files are in: ${OUTPUT_DIR}"