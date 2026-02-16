#!/bin/bash

# This code now runs as:
# ./download_files.sh /path/to/input /path/to/output
# or ./download_files.sh /path/to/input/ /path/to/output/, now that I included that helper function
# Example for this case: ./download_files_with_timeout.sh ./DownloadListFromAlien.txt /storage2/cicero/asymmetric_rapidity_V0tests_LHC25_PbPb_pass1/

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
# Tunables (speed + robustness)
########################################
PARALLEL_JOBS=2 # Lowered to 2 to keep some level of multi-jobs, while still robust.
ALIEN_TIMEOUT=5     # seconds (default is 150... too slow)
ALIEN_RETRY=1        # don’t retry broken replicas forever
# ALIEN_CP_OPTS="-x -T ${ALIEN_TIMEOUT} -retry ${ALIEN_RETRY} -t ${PARALLEL_JOBS}"
ALIEN_CP_OPTS="-T ${ALIEN_TIMEOUT} -retry ${ALIEN_RETRY} -t ${PARALLEL_JOBS}"
# PERSONAL_TIMEOUT=45

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
        # alien_cp -t 8 "alien:${FILE}" "file:${OUTFILE}"
        # alien_cp ${ALIEN_CP_OPTS} "alien:${FILE}" "file:${OUTFILE}"
        # Doing a manual timeout check to kill the process when alien_cp can't do it by hand!
        timeout 120s alien_cp -t 6 "alien:${FILE}" "file:${OUTFILE}"
        RC=$?

        if [[ $RC -eq 124 ]]; then
            echo "    120s TIMEOUT -- deleting partial file: $OUTFILE"
            rm -f "$OUTFILE"
        fi
    fi

    (( COUNTER++ ))
done < "$INPUT_FILE"