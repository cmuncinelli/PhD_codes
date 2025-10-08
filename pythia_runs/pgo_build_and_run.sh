#!/usr/bin/env bash
# =======================================================================
#  Profile-Guided Optimization (PGO) build script for PythiaGenMin
# =======================================================================
#  This script will:
#    1. Build and run a training (instrumented) version of PythiaGenMin
#    2. Build an optimized PGO+LTO version automatically afterward
#    3. Clean up all temporary profile data and intermediate binaries
#
#  Usage:
#     ./pgo_build_and_run.sh
# =======================================================================

set -euo pipefail

# === Helper to print a timestamped header ==============================
timestamp() { date +"[%H:%M:%S]"; }
header() { echo; echo "$(timestamp) === $1 ==="; }

# === Configuration =====================================================
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SRC_FILE="PythiaGenMin.cxx"
TRAIN_EXE="PythiaGenMin.profgen.exe"
FINAL_EXE="PythiaGenMin.exe"

# --- Change into the script's directory ---
cd "$SCRIPT_DIR"
echo "Working directory: $(pwd)"

if [ ! -f "$SRC_FILE" ]; then
  echo "ERROR: Source file not found at $(pwd)/$SRC_FILE"; exit 1;
fi

ARGS=(
  "/storage3/cicero/pythia_data"
  "/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg"
  "2e4"
  "50"
)

EXTRA_CFLAGS=("-I$PYTHIA8/include" "-I$FASTJET/include" $(root-config --cflags))
EXTRA_LDFLAGS=("-L$PYTHIA8/lib" "-lpythia8" "-L$FASTJET/lib" "-lfastjet" $(root-config --libs) "-fopenmp")

# =======================================================================
# 0. SETUP AND CLEANUP
# =======================================================================
header "0. Performing clean setup"
rm -f "$TRAIN_EXE" "$FINAL_EXE" ./*.gcda ./*.gcno

# =======================================================================
# 1. BUILD INSTRUMENTED BINARY
# =======================================================================
header "1. Building instrumented binary for profile generation"

g++ -std=c++17 -march=native -O2 -pipe \
    -fprofile-generate -fno-omit-frame-pointer \
    "$SRC_FILE" -o "$TRAIN_EXE" \
    "${EXTRA_CFLAGS[@]}" "${EXTRA_LDFLAGS[@]}"
echo "  -> Instrumented binary created."

# =======================================================================
# 2. RUN TRAINING WORKLOAD
# =======================================================================
header "2. Running training workload to collect profile data"

time "./$TRAIN_EXE" "${ARGS[@]}"
echo "  -> Training run completed."

# =======================================================================
# 3. VERIFY AND RENAME PROFILE DATA
# =======================================================================
header "3. Verifying and Renaming profile data"
PROFILE_FILE=$(find . -maxdepth 1 -name "*.gcda" -print -quit)

if [ -z "$PROFILE_FILE" ]; then
  echo "  [FATAL ERROR] No '.gcda' profile file was generated. Halting."
  exit 1
fi
echo "  -> Found generated profile file: $PROFILE_FILE"

# --- NEW: Rename the file to what the compiler expects ---
EXPECTED_FILE="${SRC_FILE%.*}.gcda" # This creates "PythiaGenMin.gcda"
if [ "$PROFILE_FILE" != "./$EXPECTED_FILE" ]; then
  echo "  -> Renaming '$PROFILE_FILE' to '$EXPECTED_FILE' to match compiler expectation."
  mv "$PROFILE_FILE" "$EXPECTED_FILE"
else
  echo "  -> Profile file already has the expected name."
fi

# =======================================================================
# 4. BUILD OPTIMIZED BINARY
# =======================================================================
header "4. Building optimized (PGO + LTO) binary"

g++ -std=c++17 -march=native -O2 -pipe \
    -fprofile-use -fprofile-correction \
    -flto=2 \
    "$SRC_FILE" -o "$FINAL_EXE" \
    "${EXTRA_CFLAGS[@]}" "${EXTRA_LDFLAGS[@]}"

echo "  -> Optimized binary built successfully: $FINAL_EXE"

# =======================================================================
# 5. FINAL CLEAN UP
# =======================================================================
header "5. Cleaning up intermediate files"
rm -f "$TRAIN_EXE" ./*.gcda ./*.gcno
echo "  -> Cleanup complete. Only '$FINAL_EXE' remains."