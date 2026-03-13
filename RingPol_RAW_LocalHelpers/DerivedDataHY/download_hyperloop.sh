#!/usr/bin/env bash

###############################################################################
# Script name:
#   download_hyperloop.sh
#
# Purpose:
#   Orchestrates the full download pipeline for one Hyperloop train run.
#   Given a train run ID, it:
#
#     1. Looks up the dataset name and wagon shortname in train_registry.conf.
#     2. Reads <WORK_DIR>/alien_paths.txt and uses alien_ls to discover all
#        AO2D.root and AnalysisResults.root LFNs (convert_alien_paths.sh).
#     3. Downloads all AO2D files into <WORK_DIR>/AO2Ds/.
#     4. Downloads all AnalysisResults files into a temporary staging folder.
#     5. Generates the O2 input file list (input_data_storage.txt).
#     6. Merges the AnalysisResults files into a single
#        AnalysisResults_merged.root and removes the staging folder.
#
# Usage:
#   ./download_hyperloop.sh <TRAIN_RUN_ID>
#
# Example:
#   ./download_hyperloop.sh 634083
#
# Before running:
#   - Register the train run in train_registry.conf.
#   - Run init_wagon_structure.sh to create the folder tree.
#   - Fill in <WORK_DIR>/alien_paths.txt with the comma-separated list of
#     AliEN directory paths copied from the Hyperloop interface.
#
###############################################################################

# Stop execution on first error
set -e

###############################################################################
# Absolute path to the DerivedDataHY framework directory
###############################################################################
FRAMEWORK_DIR="/home/users/cicerodm/PhD_codes/RingPol_RAW_LocalHelpers/DerivedDataHY"

###############################################################################
# Paths to helper scripts and registry (all absolute)
###############################################################################
CONVERT_ALIEN_SCRIPT="${FRAMEWORK_DIR}/convert_alien_paths.sh"
GEN_PATHS_SCRIPT="${FRAMEWORK_DIR}/gen_input_paths.sh"
DOWNLOAD_SCRIPT="${FRAMEWORK_DIR}/download_files_without_timeout.sh"
MERGE_SCRIPT="${FRAMEWORK_DIR}/merge_analysis_results.sh"
REGISTRY="${FRAMEWORK_DIR}/train_registry.conf"

###############################################################################
# Validate train run ID argument
###############################################################################
if [ -z "$1" ]; then
  echo "Usage: $0 <TRAIN_RUN_ID>"
  echo "Example: $0 634083"
  exit 1
fi

TRAIN_RUN_ID="$1"

###############################################################################
# Verify all helper scripts and registry exist
###############################################################################
for REQUIRED in "$CONVERT_ALIEN_SCRIPT" "$GEN_PATHS_SCRIPT" \
                "$DOWNLOAD_SCRIPT" "$MERGE_SCRIPT" "$REGISTRY"; do
  if [ ! -f "$REQUIRED" ]; then
    echo "Error: Required file not found: ${REQUIRED}"
    exit 1
  fi
done

chmod +x "$CONVERT_ALIEN_SCRIPT" "$GEN_PATHS_SCRIPT" \
         "$DOWNLOAD_SCRIPT" "$MERGE_SCRIPT"

###############################################################################
# Look up the train run in the registry
###############################################################################
REGISTRY_LINE=$(grep -v '^\s*#' "$REGISTRY" \
  | awk -v id="$TRAIN_RUN_ID" '$1 == id {print}')

if [ -z "$REGISTRY_LINE" ]; then
  echo "Error: Train run ID '${TRAIN_RUN_ID}' not found in ${REGISTRY}"
  echo "Please add an entry to train_registry.conf before downloading."
  exit 1
fi

DATASET_NAME=$(   echo "$REGISTRY_LINE" | awk '{print $2}')
WAGON_SHORTNAME=$(echo "$REGISTRY_LINE" | awk '{print $3}')
STORAGE_BASE=$(   echo "$REGISTRY_LINE" | awk '{print $4}')

WORK_DIR="${STORAGE_BASE}/${DATASET_NAME}/${WAGON_SHORTNAME}"

echo "============================================================"
echo "  Train run ID   : ${TRAIN_RUN_ID}"
echo "  Dataset        : ${DATASET_NAME}"
echo "  Wagon          : ${WAGON_SHORTNAME}"
echo "  Working dir    : ${WORK_DIR}"
echo "============================================================"

###############################################################################
# Sanity checks
###############################################################################
if [ ! -f "${WORK_DIR}/alien_paths.txt" ]; then
  echo "Error: alien_paths.txt not found in: ${WORK_DIR}"
  echo "Fill it with the AliEN directory paths from the Hyperloop interface."
  exit 1
fi

if ! grep -qE '[^[:space:]]' "${WORK_DIR}/alien_paths.txt"; then
  echo "Error: alien_paths.txt exists but appears to be empty."
  exit 1
fi

###############################################################################
# Step 1: Discover LFNs via alien_ls
###############################################################################
echo ""
echo "[Step 1] Discovering file paths via alien_ls..."
"$CONVERT_ALIEN_SCRIPT" "${WORK_DIR}"

###############################################################################
# Step 2: Download AO2D files into AO2Ds/
###############################################################################
echo ""
echo "[Step 2] Downloading AO2D files..."
"$DOWNLOAD_SCRIPT" \
  "${WORK_DIR}/DownloadListFromAlien.txt" \
  "${WORK_DIR}/AO2Ds" \
  "AO2D"

###############################################################################
# Step 3: Download AnalysisResults files into temp_analysis_results/
###############################################################################
echo ""
echo "[Step 3] Downloading AnalysisResults files..."
"$DOWNLOAD_SCRIPT" \
  "${WORK_DIR}/DownloadListAnalysisResults.txt" \
  "${WORK_DIR}/temp_analysis_results" \
  "AnalysisResults"

###############################################################################
# Step 4: Generate O2 input file list
###############################################################################
NUMBER_OF_AODs=$(wc -l < "${WORK_DIR}/DownloadListFromAlien.txt")
echo ""
echo "[Step 4] Generating O2 input list (${NUMBER_OF_AODs} files)..."
"$GEN_PATHS_SCRIPT" "${WORK_DIR}" "${NUMBER_OF_AODs}"

###############################################################################
# Step 5: Merge AnalysisResults and clean up staging folder
###############################################################################
echo ""
echo "[Step 5] Merging AnalysisResults files..."
"$MERGE_SCRIPT" "${WORK_DIR}"

echo ""
echo "============================================================"
echo "  Done."
echo "  AO2D files     : ${WORK_DIR}/AO2Ds/"
echo "  O2 input list  : ${WORK_DIR}/input_data_storage.txt"
echo "  Merged results : ${WORK_DIR}/AnalysisResults_merged.root"
echo "============================================================"
