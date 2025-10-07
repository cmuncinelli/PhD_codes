#!/bin/bash
# # extract_top_functions_fast.sh
# # Usage: ./extract_top_functions_fast.sh perf.data output.txt

# PERF_DATA=${1:-perf.data}          # Default to perf.data if not provided
# OUTPUT_FILE=${2:-top_functions.txt}  # Default output file

# PROGRAM_NAME="PythiaGenMin.exe"

# echo "Processing perf data for program: $PROGRAM_NAME"
# echo "Output will be saved to: $OUTPUT_FILE"
# echo ""

# # Method 1: Using perf report in fast mode (no full symbol resolution)
# echo "=== Top functions using perf report (stdio, fast) ===" > "$OUTPUT_FILE"
# perf report -i "$PERF_DATA" --stdio --comm "$PROGRAM_NAME" --sort overhead -n | head -n 40 >> "$OUTPUT_FILE"

# echo "" >> "$OUTPUT_FILE"
# echo "=== Top functions using perf script + awk (efficient) ===" >> "$OUTPUT_FILE"

# # Method 2: Efficient perf script processing
# # Extract only function names, aggregate with awk, and sort
# perf script -i "$PERF_DATA" --comm "$PROGRAM_NAME" --ns \
#     | awk '{counts[$6]++} END {for (f in counts) print counts[f], f}' \
#     | sort -nr | head -20 >> "$OUTPUT_FILE"

# echo ""
# echo "Done! Top functions exported to $OUTPUT_FILE"






# Additional parallelization:
# extract_top_functions_parallel.sh
# Usage: ./extract_top_functions_parallel.sh perf.data output.txt [num_threads] [func_col]
#
# - perf.data: path to perf data (default: perf.data)
# - output.txt: output file (default: top_functions.txt)
# - num_threads: optional, default = nproc
# - func_col: optional, awk column that contains function name, default = 6
#
# Notes:
#  - If GNU parallel is present, this uses `parallel --pipe` (preferred, streaming).
#  - If not, it falls back to splitting the perf script output into chunk files and
#    processing them with xargs -P.
#  - Final aggregation is done with awk and sort (no huge in-memory sorts).
#  - Choose NUM_THREADS sensibly (nproc is a safe default).

set -euo pipefail

PERF_DATA="${1:-perf.data}"
OUTPUT_FILE="${2:-top_functions.txt}"
NUM_THREADS="${3:-128}"
FUNC_COL="${4:-6}"        # column where perf script places the function name
BLOCK_SIZE="${BLOCK_SIZE:-8M}"  # block size for parallel --pipe (8M default)

# Ensure output folder:
OUT_DIR="profiling"
OUTPUT_FILE="$OUT_DIR/$OUTPUT_FILE"

# Safety checks
if [[ ! -f "$PERF_DATA" ]]; then
  echo "ERROR: perf data file not found: $PERF_DATA" >&2
  exit 2
fi

# Make NUM_THREADS an integer and cap it to nproc
NPROC="$(nproc --all)"
if ! [[ "$NUM_THREADS" =~ ^[0-9]+$ ]]; then
  NUM_THREADS="$NPROC"
fi
if (( NUM_THREADS > NPROC )); then
  echo "Warning: requested threads ($NUM_THREADS) > available cores ($NPROC). Using $NPROC."
  NUM_THREADS="$NPROC"
fi

echo "Perf file: $PERF_DATA"
echo "Program filter (comm): PythiaGenMin.exe"
echo "Output file: $OUTPUT_FILE"
echo "Using up to $NUM_THREADS worker(s)"
echo "Using func_col=$FUNC_COL"
echo ""

# Header
{
  echo "Top functions for program: PythiaGenMin.exe"
  echo "Perf data: $PERF_DATA"
  echo "Generated: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
  echo ""
  echo "=== Aggregation (top 20) ==="
} > "$OUTPUT_FILE"

# Helper: merge partial counts (tab-separated "count<TAB>func")
merge_and_output() {
  # Input: stream of lines "count<TAB>func"
  # Output: top 20 lines "count<TAB>func" sorted by count desc
  awk -F'\t' '{ key=$2; if (key=="") next; counts[key]+=$1 } END { for (k in counts) printf("%d\t%s\n", counts[k], k) }' \
    | sort -nr -k1,1 | head -n 20
}

# Preferred method: GNU parallel --pipe (no temp chunk files)
if command -v parallel >/dev/null 2>&1; then
  echo "Using GNU parallel (streaming) mode..."
  # perf script -> parallel --pipe (splits into blocks) -> awk per-block -> merge
  # Use --block to control block size (default 8M). The number of concurrent jobs is -j.
  # We print count<TAB>func from each worker and then merge.
  perf script -i "$PERF_DATA" --comm "PythiaGenMin.exe" --ns 2>/dev/null \
    | parallel --jobs "$NUM_THREADS" --pipe --block "$BLOCK_SIZE" \
        awk -v col="$FUNC_COL" '{
            # crude extraction: prefer column "col", fall back to "UNKNOWN"
            if (NF >= col) f = $col; else f = "UNKNOWN";
            counts[f]++
        }
        END {
            for (k in counts) printf("%d\t%s\n", counts[k], k)
        }' \
    | merge_and_output >> "$OUTPUT_FILE"

  echo "Done (parallel mode). Results appended to $OUTPUT_FILE"
  exit 0
fi

# If we get here, GNU parallel is not available. Try split + xargs fallback.
echo "GNU parallel not found. Falling back to split + xargs mode."

# Create tmp dir and ensure cleanup
TMPDIR="$(mktemp -d -t extract_perf_chunks.XXXXXX)"
cleanup() {
  rm -rf "$TMPDIR"
}
trap cleanup EXIT

# Decide a chunk byte size for `split`. Use perf.data size / NUM_THREADS (min 8M).
PERF_BYTES=$(stat -c%s "$PERF_DATA")
MIN_BYTES=$((8 * 1024 * 1024))  # 8MB
CHUNK_BYTES=$(( (PERF_BYTES / NUM_THREADS) + 1 ))
if (( CHUNK_BYTES < MIN_BYTES )); then
  CHUNK_BYTES=$MIN_BYTES
fi

echo "Creating chunk files in $TMPDIR (approx chunk bytes: $CHUNK_BYTES)..."
# Stream perf script output into split to create chunk files (no single enormous temp file)
# split -b accepts bytes; prefix is $TMPDIR/chunk_
perf script -i "$PERF_DATA" --comm "PythiaGenMin.exe" --ns 2>/dev/null \
  | split -b "$CHUNK_BYTES" - "$TMPDIR/chunk_"

CHUNK_COUNT=$(ls "$TMPDIR"/chunk_* 2>/dev/null | wc -l || echo 0)
if (( CHUNK_COUNT == 0 )); then
  echo "ERROR: no chunk files were produced; falling back to single-threaded mode." >&2
  # fall through to single-threaded below
else
  echo "Processing $CHUNK_COUNT chunk(s) with up to $NUM_THREADS parallel workers..."

  # Create a helper chunk processor script so we avoid quoting hell
  cat > "$TMPDIR/process_chunk.sh" <<'SH'
#!/usr/bin/env bash
COL="$1"
FILE="$2"
awk -v col="$COL" '{
  if (NF >= col) f = $col; else f = "UNKNOWN";
  counts[f]++
}
END {
  for (k in counts) printf("%d\t%s\n", counts[k], k)
}' "$FILE"
SH
  chmod +x "$TMPDIR/process_chunk.sh"

  # Use xargs to run the chunk processor in parallel and aggregate
  ls "$TMPDIR"/chunk_* \
    | xargs -n 1 -P "$NUM_THREADS" -I {} "$TMPDIR/process_chunk.sh" "$FUNC_COL" {} \
    | merge_and_output >> "$OUTPUT_FILE"

  echo "Done (split+xargs mode). Results appended to $OUTPUT_FILE"
  exit 0
fi

# Final fallback: single-threaded streaming awk (works but may be slower)
echo "Falling back to single-threaded aggregation (awk). This may be slow for very large perf.data."
perf script -i "$PERF_DATA" --comm "PythiaGenMin.exe" --ns 2>/dev/null \
  | awk -v col="$FUNC_COL" '{
      if (NF >= col) f = $col; else f = "UNKNOWN";
      counts[f]++
    }
    END {
      for (k in counts) printf("%d\t%s\n", counts[k], k)
    }' \
  | sort -nr -k1,1 | head -n 20 >> "$OUTPUT_FILE"

echo "Done (single-threaded). Results appended to $OUTPUT_FILE"
exit 0
