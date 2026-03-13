#!/usr/bin/env bash

###############################################################################
# Script name:
#   merge_analysis_results.sh
#
# Purpose:
#   Merges all individual AnalysisResults_N.root files downloaded from AliEN
#   into a single file using ROOT's hadd, then deletes the temporary folder
#   they were staged in.
#
#   Input  : <WORK_DIR>/temp_analysis_results/AnalysisResults_*.root
#   Output : <WORK_DIR>/AnalysisResults_merged.root
#
#   After a successful merge the temp_analysis_results/ directory is removed.
#   If hadd fails, the temp directory is preserved so you can inspect or
#   retry without re-downloading.
#
# Arguments:
#   $1 : WORK_DIR -- absolute path to the per-wagon working directory
#
# Example:
#   ./merge_analysis_results.sh \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2_small/ITSandTPC
#
# Dependencies:
#   hadd -- provided by the ROOT / O2 environment. Make sure the O2
#           environment is loaded before calling this script.
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

TEMP_DIR="${WORK_DIR}/temp_analysis_results"
OUTPUT_FILE="${WORK_DIR}/AnalysisResults_merged.root"

###############################################################################
# Check that the temp directory exists and contains files to merge
###############################################################################
if [ ! -d "$TEMP_DIR" ]; then
  echo "Error: Staging directory not found: ${TEMP_DIR}"
  echo "Run the AnalysisResults download step first."
  exit 1
fi

# Build a sorted list of files to merge (version sort handles _10 > _9 correctly)
shopt -s nullglob
INPUT_FILES=("${TEMP_DIR}"/AnalysisResults_*.root)
shopt -u nullglob

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
  echo "Error: No AnalysisResults_*.root files found in: ${TEMP_DIR}"
  exit 1
fi

echo "============================================================"
echo "  Merging ${#INPUT_FILES[@]} AnalysisResults file(s)"
echo "  Input  : ${TEMP_DIR}/AnalysisResults_*.root"
echo "  Output : ${OUTPUT_FILE}"
echo "============================================================"

###############################################################################
# Run hadd
#
# hadd arguments:
#   -f : overwrite output file if it already exists
#
# The input files are passed as a sorted glob expansion (shell-sorted by
# the array construction above, which uses version sort via nullglob +
# the natural order of the glob).
###############################################################################
hadd -f "$OUTPUT_FILE" "${INPUT_FILES[@]}"

HADD_EXIT=$?

###############################################################################
# Stage-out: remove temp directory only on success
###############################################################################
if [ $HADD_EXIT -eq 0 ]; then
  echo ""
  echo "  [OK] hadd succeeded."
  echo "  Removing staging directory: ${TEMP_DIR}"
  rm -rf "$TEMP_DIR"
  echo ""
  echo "SUCCESS."
  echo "  Merged output : ${OUTPUT_FILE}"
else
  echo ""
  echo "  [FAIL] hadd exited with code ${HADD_EXIT}."
  echo "  Staging directory preserved for inspection: ${TEMP_DIR}"
  exit $HADD_EXIT
fi
