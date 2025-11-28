#!/bin/bash

########################################
# Resolve script directory
########################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

########################################
# Helper: remove trailing slashes
########################################
strip_trailing_slashes() {
    local path="$1"
    # Remove all trailing slashes using parameter expansion
    path="${path%/}"
    echo "$path"
}

########################################
# Input arguments or defaults
########################################
RAW_INPUT_DIR="${1:-$SCRIPT_DIR}"
RAW_OUTPUT_DIR="${2:-$SCRIPT_DIR}"

# Normalize paths: strip trailing slashes
INPUT_DIR=$(strip_trailing_slashes "$RAW_INPUT_DIR")
OUTPUT_DIR=$(strip_trailing_slashes "$RAW_OUTPUT_DIR")

INPUT_FILE="${INPUT_DIR}/DownloadListFromAlien.txt"

########################################
# Basic sanity checks
########################################
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

########################################
# Main loop
########################################

COUNTER=1
MAX=$(wc -l < "$INPUT_FILE")

while IFS= read -r FILE; do
    echo "[FILE #${COUNTER}/${MAX}: ${FILE}]"

    OUTFILE="${OUTPUT_DIR}/AO2D_${COUNTER}.root"

    if [[ -f "$OUTFILE" ]]; then
        echo "File AO2D_${COUNTER}.root already exists! Skipping."
    else
        echo "Copying: alien:${FILE}  ->  $OUTFILE"
        alien_cp -t 8 "alien:${FILE}" "file:${OUTFILE}"
    fi

    (( COUNTER++ ))
done < "$INPUT_FILE"

# This code now runs as:
# ./download_files.sh /path/to/input /path/to/output
# or ./download_files.sh /path/to/input/ /path/to/output/, now that I included that helper function