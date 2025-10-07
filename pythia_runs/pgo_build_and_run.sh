#!/usr/bin/env bash
# =======================================================================
#  Profile-Guided Optimization (PGO) build script for PythiaGenMin
# =======================================================================
#  Usage:
#     ./pgo_build_and_run.sh train   # build with -fprofile-generate and run workload
#     ./pgo_build_and_run.sh final   # build optimized binary using collected profiles
#
#  Notes:
#   - Training run produces profiling data in ./pgodata/
#   - The final binary will use -O2, LTO, and PGO data for optimization
# =======================================================================

set -e  # exit immediately on any error
set -u  # treat unset variables as errors

# === Configuration =====================================================
SRC="PythiaGenMin.cxx"
TRAIN_DIR="./pgodata"

# === Program arguments (safely stored in array) ========================
ARGS=(
  "/storage3/cicero/pythia_data"
  "/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg"
  "1e5"
  "324"
)

# === Common includes and libraries =====================================
EXTRA_CFLAGS=(
  "-I$PYTHIA8/include"
  "-I$FASTJET/include"
  $(root-config --cflags)
)

EXTRA_LDFLAGS=(
  "-L$PYTHIA8/lib" "-lpythia8"
  "-L$FASTJET/lib" "-lfastjet"
  $(root-config --libs)
  "-fopenmp"
)

# =======================================================================
# Build for profile generation (instrumented binary)
# =======================================================================
if [[ "${1:-}" == "train" ]]; then
  echo "=== Building instrumented (profile generation) binary ==="

  g++ -std=c++17 -march=native -O2 -pipe \
      -fprofile-generate=$TRAIN_DIR -fno-omit-frame-pointer \
      "$SRC" -o PythiaGenMin.profgen.exe \
      "${EXTRA_CFLAGS[@]}" "${EXTRA_LDFLAGS[@]}"

  echo
  echo "=== Running training workload (this will generate profile data) ==="
  echo "Command line:"
  echo "  ./PythiaGenMin.profgen.exe \"${ARGS[@]}\""
  echo

  ./PythiaGenMin.profgen.exe "${ARGS[@]}"

  echo
  echo "=== Training run complete ==="
  echo "Profile data stored in: $TRAIN_DIR"
  echo "Next step: ./pgo_build_and_run.sh final"
  exit 0
fi

# =======================================================================
# Build optimized binary using collected profile data
# =======================================================================
if [[ "${1:-}" == "final" ]]; then
  if [[ ! -d "$TRAIN_DIR" ]]; then
    echo "Error: Profile data directory '$TRAIN_DIR' not found. Run 'train' first."
    exit 1
  fi

  echo "=== Building optimized (PGO + LTO) binary ==="

  g++ -std=c++17 -march=native -O2 -pipe \
      -fprofile-use=$TRAIN_DIR -fprofile-correction \
      -flto=2 \
      "$SRC" -o PythiaGenMin.exe \
      "${EXTRA_CFLAGS[@]}" "${EXTRA_LDFLAGS[@]}"

  echo
  echo "=== Optimized binary built: ./PythiaGenMin.exe ==="
  echo "To test:"
  echo "  time ./PythiaGenMin.exe \"${ARGS[@]}\""
  exit 0
fi

# =======================================================================
# Help message
# =======================================================================
echo "Usage:"
echo "  $0 train    # build and run the instrumented version (collect profiles)"
echo "  $0 final    # build optimized binary using the profile data"
exit 1
