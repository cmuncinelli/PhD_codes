#!/usr/bin/env bash

###############################################################################
# Script name:
#   convert_alien_paths.sh
#
# Purpose:
#   Reads a list of AliEN directory paths from:
#
#     <WORK_DIR>/alien_paths.txt
#
#   Each path points to a Hyperloop output folder (e.g. .../hy_XXXXXXX/AOD).
#   The script uses alien_ls to enumerate the one-level-deep subdirectories
#   inside each path (e.g. 001, 002, ...) and constructs the full LFN for
#   the two files expected in each subdirectory:
#
#     AO2D.root            -> appended to DownloadListFromAlien.txt
#     AnalysisResults.root -> appended to DownloadListAnalysisResults.txt
#
#   Both output files are written to <WORK_DIR> and overwrite any previous
#   version.
#
# Input format (alien_paths.txt):
#   Paths may be comma-separated on one line (paste directly from HY),
#   one per line, or any mix of both. Blank lines and surrounding whitespace
#   are ignored.
#
#   Example (single pasted line from HY):
#     /alice/.../hy_4387490/AOD,/alice/.../hy_4387489/AOD
#
#   Example (one per line):
#     /alice/.../hy_4387490/AOD
#     /alice/.../hy_4387489/AOD
#
# Arguments:
#   $1 : WORK_DIR -- absolute path to the per-wagon working directory
#
# Example:
#   ./convert_alien_paths.sh \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2_small/ITSandTPC
#
###############################################################################

# Stop on first error
set -e

###############################################################################
# Read and validate WORK_DIR
###############################################################################
WORK_DIR="$1"

if [ -z "$WORK_DIR" ]; then
  echo "Usage: $0 <WORK_DIR>"
  exit 1
fi

if [ ! -d "$WORK_DIR" ]; then
  echo "Error: WORK_DIR does not exist: ${WORK_DIR}"
  exit 1
fi

ALIEN_PATHS_FILE="${WORK_DIR}/alien_paths.txt"
AOD_LIST="${WORK_DIR}/DownloadListFromAlien.txt"
AR_LIST="${WORK_DIR}/DownloadListAnalysisResults.txt"

if [ ! -f "$ALIEN_PATHS_FILE" ]; then
  echo "Error: alien_paths.txt not found in: ${WORK_DIR}"
  echo "Fill in the file created by init_wagon_structure.sh before running."
  exit 1
fi

###############################################################################
# Parse alien_paths.txt
#
# Strategy:
#   1. Replace all commas with newlines so the paste-as-is format is handled.
#   2. Trim leading/trailing whitespace from each resulting token.
#   3. Drop blank lines.
#
# The result is a clean newline-separated list of AliEN directory paths
# stored in the PATHS variable.
###############################################################################
PATHS=$(tr ',' '\n' < "$ALIEN_PATHS_FILE" \
  | sed 's/^[[:space:]]*//' \
  | sed 's/[[:space:]]*$//' \
  | grep -v '^$')

if [ -z "$PATHS" ]; then
  echo "Error: No paths found in ${ALIEN_PATHS_FILE}."
  exit 1
fi

###############################################################################
# Clear previous output files
###############################################################################
rm -f "$AOD_LIST" "$AR_LIST"

###############################################################################
# Main loop: alien_ls each provided directory, then construct LFNs
#
# Expected structure (one level deep, confirmed by user):
#   <ALIEN_DIR>/001/AO2D.root
#   <ALIEN_DIR>/001/AnalysisResults.root
#   <ALIEN_DIR>/002/AO2D.root
#   ...
###############################################################################
TOTAL_PATHS=0
TOTAL_FILES=0

while IFS= read -r ALIEN_DIR || [ -n "$ALIEN_DIR" ]; do

  # Skip blank lines that may survive the earlier filtering
  [ -z "$ALIEN_DIR" ] && continue

  (( TOTAL_PATHS += 1 ))

  # alien_ls returns one entry name per line (just the basename, not full path).
  # The "|| true" prevents set -e from aborting the script if alien_ls exits
  # non-zero, which some O2 environments do even on a successful listing.
  SUBDIRS=$(alien_ls "${ALIEN_DIR}" 2>&1) || true

  if [ -z "$SUBDIRS" ]; then
    echo "  Warning: alien_ls returned nothing for ${ALIEN_DIR} -- skipping."
    continue
  fi

  # If alien_ls returned an error message instead of a directory listing,
  # print it so the user can see what went wrong.
  if echo "$SUBDIRS" | grep -qi "error\|fail\|denied\|no such"; then
    echo "  Warning: alien_ls may have returned an error for ${ALIEN_DIR}:"
    echo "    ${SUBDIRS}"
  fi

  while IFS= read -r SUBDIR || [ -n "$SUBDIR" ]; do

    [ -z "$SUBDIR" ] && continue

    # alien_ls sometimes returns non-directory entries such as
    # "aod_collection.xml". Skip anything that is not a pure numeric
    # subdirectory name (with or without a trailing slash).
    CLEAN="${SUBDIR%/}"
    if ! [[ "$CLEAN" =~ ^[0-9]+$ ]]; then
      continue
    fi

    # alien_ls may return entries with a trailing slash (e.g. "001/").
    # Strip it so the concatenated path does not contain a double slash.
    SUBDIR="${CLEAN}"

    echo "${ALIEN_DIR}/${SUBDIR}/AO2D.root"            >> "$AOD_LIST"
    echo "${ALIEN_DIR}/${SUBDIR}/AnalysisResults.root" >> "$AR_LIST"
    (( TOTAL_FILES += 1 ))

  done <<< "$SUBDIRS"

done <<< "$PATHS"

###############################################################################
# Summary
###############################################################################
echo ""
echo "Done."
echo "  Processed    : ${TOTAL_PATHS} AliEN director(ies)"
echo "  File pairs   : ${TOTAL_FILES} (AO2D + AnalysisResults each)"
echo "  AOD list     : ${AOD_LIST}"
echo "  AR list      : ${AR_LIST}"