#!/bin/bash

# This command merges the previous list I made by hand with the new list, without changing the order of the old list.
# The download script needs to see the same order as before, so we need to preserve the old list's order, or be forced
# to redownload everything (not happening!!!)

OLD="${1:-DownloadListFromAlien_saved.txt}"
NEW="${2:-DownloadListFromAlien.txt}"
OUT="${3:-DownloadListFromAlien.txt}"

# Basic sanity checks
if [[ ! -f "$OLD" ]]; then
    echo "ERROR: Old file not found: $OLD"
    exit 1
fi

if [[ ! -f "$NEW" ]]; then
    echo "ERROR: New file not found: $NEW"
    exit 1
fi

# Temp file to store the merged list
TMP=$(mktemp)

# Step 1 — Add all old entries in original order
cat "$OLD" >> "$TMP"

# Ensure TMP ends with a newline so appended lines don't join the last line
if [[ -s "$TMP" && "$(tail -c1 "$TMP")" != $'\n' ]]; then
    printf '\n' >> "$TMP"
fi

# Step 2 — Append new entries only if they were not previously present
grep -vxFf "$OLD" "$NEW" >> "$TMP"

# Step 3 — Move merged result to output file
mv "$TMP" "$OUT"

echo "Merged list written to $OUT"
