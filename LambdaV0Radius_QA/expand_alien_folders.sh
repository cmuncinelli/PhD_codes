#!/bin/bash

# Create a list of FoldersToExpand from Hyperloop's "Copy all output directories" list for a given train run,
# then just run this code to unfold everything and then proceed to the "download_files.sh" script and the
# "input_data_storage.txt" script.

INPUT="${1:-FoldersToExpand.txt}"
OUTPUT="${2:-DownloadListFromAlien.txt}"

# Check input file
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

# Count total lines (for progress)
TOTAL=$(grep -cve '^\s*$' "$INPUT")
COUNT=1

# Empty output file
> "$OUTPUT"

while IFS= read -r DIR; do
    # Skip empty lines
    [[ -z "$DIR" ]] && continue

    # # Skip entries WITHOUT trailing "/AOD"
    # if [[ ! "$DIR" =~ /AOD$ ]]; then
    #     echo "Skipping (not an AOD folder): $DIR"
    #     (( COUNT++ ))
    #     continue
    # fi
    # No need to do these skips! There are valid files like:
    ### /alice/cern.ch/user/a/alihyperloop/jobs/0308/hy_3088070/0001/AO2D.root
    # in the repo! No need to skip those folders!


    PERCENT=$(( 100 * COUNT / TOTAL ))
    echo -e "Expanding: $DIR\t [${COUNT}/${TOTAL}] (${PERCENT}%)"

    # List subdirs like 001, 002, 010, ...
    SUBDIRS=$(alien_ls "$DIR" | grep -E '^[0-9]+/?$')

    for S in $SUBDIRS; do
        S="${S%/}"
        FILE="${DIR}/${S}/AO2D.root"
        echo "$FILE" >> "$OUTPUT"
    done

    (( COUNT++ ))

done < "$INPUT"
