#!/usr/bin/env bash

###############################################################################
# Script name:
#   init_wagon_structure.sh
#
# Purpose:
#   Reads a registry configuration file and creates the full per-wagon directory
#   structure. Safe to re-run at any time: existing directories and files are 
#   never overwritten or deleted.
#
#   For each entry in the registry, the following tree is created:
#
#     <STORAGE_BASE>/<DATASET_NAME>/<WAGON_SHORTNAME>/
#       aod_collections/           <- place MonALISA XML files here manually
#       AO2Ds/                     <- AO2D_*.root files land here
#       results_consumer/          <- output of runDerivedDataConsumer_HY.sh
#       results_consumer/logs/
#       results_consumer/used_configs/
#
# Usage:
#   ./init_wagon_structure.sh [REGISTRY_FILE]
#
#   REGISTRY_FILE is optional. If omitted, the script defaults to using the 
#   legacy 'train_registry.conf' located in the same directory as this script.
#
# Example (Legacy / Hyperloop workflow):
#   ./init_wagon_structure.sh
#
# Example (Local offline workflow):
#   ./init_wagon_structure.sh /home/users/cicerodm/RingPol/LHC25ae_pass2_local/local_registry_auto.conf
#
###############################################################################

# Stop execution on first error
set -e

###############################################################################
# Configure paths and read arguments
###############################################################################
# Absolute path to this script's directory (the DerivedDataHY repo folder)
FRAMEWORK_DIR="/home/users/cicerodm/PhD_codes/RingPol_RAW_LocalHelpers/DerivedDataHY"
DEFAULT_REGISTRY="${FRAMEWORK_DIR}/train_registry.conf"

# Use the first argument if provided, otherwise fall back to the legacy default
REGISTRY="${1:-$DEFAULT_REGISTRY}"

if [ ! -f "$REGISTRY" ]; then
  echo "Error: Registry not found at: ${REGISTRY}"
  exit 1
fi

echo "========================================================"
echo "  Initializing Wagon Structures (Consumer)"
echo "  Registry: ${REGISTRY}"
echo "========================================================"
echo ""

CREATED=0
SKIPPED=0

###############################################################################
# Process each non-comment, non-empty line in the registry
###############################################################################
while IFS= read -r LINE || [ -n "$LINE" ]; do

  # Skip comment lines and blank lines
  [[ "$LINE" =~ ^\s*# ]] && continue
  [[ -z "${LINE// }" ]]  && continue

  # Parse the four fields expected in the registry format
  # Format: TRAIN_RUN_ID  DATASET_NAME  WAGON_SHORTNAME  STORAGE_BASE
  TRAIN_RUN_ID=$( echo "$LINE" | awk '{print $1}')
  DATASET_NAME=$(  echo "$LINE" | awk '{print $2}')
  WAGON_SHORTNAME=$(echo "$LINE" | awk '{print $3}')
  STORAGE_BASE=$(  echo "$LINE" | awk '{print $4}')

  WORK_DIR="${STORAGE_BASE}/${DATASET_NAME}/${WAGON_SHORTNAME}"

  echo "  Train run ${TRAIN_RUN_ID}  ->  ${WORK_DIR}"

  # Create all required subdirectories (mkdir -p is a no-op if they exist)
  # This guarantees the consumer script will have a place to write its logs
  # and results, regardless of whether the AO2Ds were downloaded or produced locally.
  mkdir -p "${WORK_DIR}/aod_collections"
  mkdir -p "${WORK_DIR}/AO2Ds"
  mkdir -p "${WORK_DIR}/results_consumer/logs"
  mkdir -p "${WORK_DIR}/results_consumer/used_configs"

  # Create a blank alien_paths.txt if it does not yet exist.
  # Paste the comma-separated AliEN directory paths from HY into this
  # file before running download_hyperloop.sh.
  if [ ! -f "${WORK_DIR}/alien_paths.txt" ]; then
    touch "${WORK_DIR}/alien_paths.txt"
  fi

  # Count newly created vs already existing using a sentinel file
  SENTINEL="${WORK_DIR}/.initialized_consumer"
  if [ ! -f "$SENTINEL" ]; then
    # Record the registry details and creation timestamp in the sentinel
    echo "train_run_id=${TRAIN_RUN_ID}" >  "$SENTINEL"
    echo "dataset=${DATASET_NAME}"      >> "$SENTINEL"
    echo "wagon=${WAGON_SHORTNAME}"     >> "$SENTINEL"
    echo "created=$(date -Iseconds)"    >> "$SENTINEL"
    echo "    Created consumer structure."
    (( ++CREATED ))
  else
    echo "    Already initialized -- skipped."
    (( ++SKIPPED ))
  fi

done < "$REGISTRY"

echo ""
echo "Done.  Created: ${CREATED}  |  Already existed: ${SKIPPED}"
echo "========================================================"