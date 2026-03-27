#!/usr/bin/env bash
###############################################################################
# Script name:
#   init_localWagon_structure.sh
#
# Purpose:
#   Pre-scaffolds the local directory structure for every wagon you plan to 
#   produce. Instead of a manual registry, it automatically infers the wagons 
#   by scanning your producer configs directory.
#
#   For each dpl-config-*.json found, it creates:
#     <OUTPUT_DATASET_DIR>/<WAGON_SHORTNAME>/
#       AO2Ds/                         <- Local AODs will land here
#       results_producer/logs/         <- O2 pipeline batch logs
#       results_producer/used_configs/ <- Snapshots of the JSON configs
#
# Usage:
#   ./init_localWagon_structure.sh [OUTPUT_DATASET_DIR] [PRODUCER_CONFIGS_DIR]
#
###############################################################################

# Stop execution on first error
set -e

# ==============================================================================
# DEFAULTS & ARGUMENTS
# ==============================================================================
DEFAULT_DATASET_DIR="/home/users/cicerodm/RingPol/LHC25ae_pass2_local"
DEFAULT_CONFIGS_DIR="/home/users/cicerodm/RingPol/producer_configs"

OUTPUT_DATASET_DIR="${1:-$DEFAULT_DATASET_DIR}"
PRODUCER_CONFIGS_DIR="${2:-$DEFAULT_CONFIGS_DIR}"

if [ ! -d "$PRODUCER_CONFIGS_DIR" ]; then
  echo "Error: Configs directory not found: ${PRODUCER_CONFIGS_DIR}"
  exit 1
fi

echo "========================================================"
echo "  Initializing Local Wagon Structures"
echo "  Dataset Dir : ${OUTPUT_DATASET_DIR}"
echo "  Configs Dir : ${PRODUCER_CONFIGS_DIR}"
echo "========================================================"
echo ""

CREATED=0
SKIPPED=0

shopt -s nullglob
CONFIG_FILES=("${PRODUCER_CONFIGS_DIR}"/dpl-config-*.json)
shopt -u nullglob

if [ ${#CONFIG_FILES[@]} -eq 0 ]; then
  echo "No dpl-config-*.json files found in ${PRODUCER_CONFIGS_DIR}. Nothing to do."
  exit 0
fi

# ==============================================================================
# PROCESS CONFIGS
# ==============================================================================
for CONFIG_FILE in "${CONFIG_FILES[@]}"; do
  
  # Extract WAGON_SHORTNAME (e.g., "dpl-config-ITSandTPC.json" -> "ITSandTPC")
  CONFIG_BASENAME=$(basename "$CONFIG_FILE" .json)
  WAGON_SHORTNAME="${CONFIG_BASENAME#dpl-config-}"
  
  WORK_DIR="${OUTPUT_DATASET_DIR}/${WAGON_SHORTNAME}"

  echo "  Scaffolding Wagon: ${WAGON_SHORTNAME} -> ${WORK_DIR}"

  # Create all required subdirectories
  mkdir -p "${WORK_DIR}/AO2Ds"
  mkdir -p "${WORK_DIR}/results_producer/logs"
  mkdir -p "${WORK_DIR}/results_producer/used_configs"

  # Count newly created vs already existing using the sentinel file
  SENTINEL="${WORK_DIR}/.initialized_producer"
  if [ ! -f "$SENTINEL" ]; then
    echo "wagon=${WAGON_SHORTNAME}"   > "$SENTINEL"
    echo "created=$(date -Iseconds)" >> "$SENTINEL"
    echo "    Created."
    (( ++CREATED ))
  else
    echo "    Already initialized -- skipped."
    (( ++SKIPPED ))
  fi

done

echo ""
echo "Done.  Created: ${CREATED}  |  Already existed: ${SKIPPED}"
echo "========================================================"