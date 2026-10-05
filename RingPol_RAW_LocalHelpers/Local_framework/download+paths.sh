#!/usr/bin/env bash

# Stop execution on first error
set -e

# Make sure all scripts are executable
chmod +x *.sh

# Base directory where AO2D files will be downloaded
STORAGE_PATH="/storage2/cicero/RingPol/LHC25an_pass1_small_partialRawDownload"

# Step 1:
# Convert XML collections into a plain list of LFNs
/home/users/cicerodm/PhD_codes/RingPol_RAW_LocalHelpers/Local_framework/convert_list_of_paths.sh

# Step 2:
# Count how many AO2D entries were produced
# This assumes one path per line in DownloadListFromAlien.txt
NUMBER_OF_AODs=$(wc -l < DownloadListFromAlien.txt)

echo "Number of AO2D files detected: ${NUMBER_OF_AODs}"

# Step 3:
# Generate the local input file list for O2
/home/users/cicerodm/PhD_codes/RingPol_RAW_LocalHelpers/Local_framework/gen_input_paths.sh "$STORAGE_PATH" "$NUMBER_OF_AODs"

# Step 4:
# Download files from ALIEN to local storage
# /home/users/cicerodm/PhD_codes/RingPol_RAW_LocalHelpers/download_files_with_timeout.sh ./ "$STORAGE_PATH"
    # Alternative code, without timeout:
    ## An alternative setting, without the timeout, to deal with the RAW AODs being way too big (3GB each, instead of 300MB)
/home/users/cicerodm/PhD_codes/RingPol_RAW_LocalHelpers/Local_framework/download_files_without_timeout.sh ./ "$STORAGE_PATH"