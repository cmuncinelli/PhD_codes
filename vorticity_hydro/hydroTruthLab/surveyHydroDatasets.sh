#!/usr/bin/env bash

###############################################################################
# Script name:
#   surveyHydroDatasets.sh
#
# Purpose:
#   Read-only survey of Vitor's hydro folders, to fill the datasets registry of runHydroTruthLab.sh from facts. For every dataset folder:
#
#     - its top: subfolders and files, and the parameter files (.dat) with their first two lines and the number of fields of their rows
#     - every sample, i.e. a folder whose numbered subfolders are events, or a single-event folder (one with a MUSIC input of its own):
#       the event IDs as ranges, how many events hold the iSS file, the grid of the first one, and the MUSIC inputs (initial condition and
#       bullet block), each key either constant, or with its number of distinct values and their range
#     - the pairing of every no_bullet-like sample with every other sample of the same folder: the same initial condition, event by event.
#       Initial conditions are compared through the last three components of their path (<folder>/<event>/<file>), since some folders were
#       moved between machines: a pair that differs only in the prefix is counted as the same and reported as such
#     - with --pair, the same pairing for any two samples, also across dataset folders
#
#   Repeated information is summarized, never listed: event IDs as ranges, values as "constant" or ranges, a grid identical to an earlier
#   one as "same as", long lists cut with "(+N more)". Nothing is written outside a temporary folder, which is removed at exit.
#
# Usage:
#   ./surveyHydroDatasets.sh [OPTIONS] [DATASET ...]
#
# Options:
#   -h, --help                 Show this help message.
#   -b, --base DIR             Folder holding the datasets (default: /storage1/vribeiro). DATASET names are relative to it.
#   -i, --iss NAME             iSS file looked for in each event (default: Smu_dpTdphidy_Thermal_rapidity_3122.dat).
#   -p, --pair SAMPLE NOBULLET Surveys both samples (paths relative to the base, or absolute) and compares their initial conditions.
#                              Repeatable.
#
#   With no DATASET and no --pair, every folder under the base is surveyed.
#
# Examples:
#   ./surveyHydroDatasets.sh > surveyHydroDatasets.log
#   ./surveyHydroDatasets.sh lumpy_events_40_50 lumpy_events_60_70 smooth_events
#   ./surveyHydroDatasets.sh --pair position_lumpy_events_40_50/R+2 lumpy_events_40_50/no_bullet
#
###############################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE="/storage1/vribeiro"
ISS_NAME="Smu_dpTdphidy_Thermal_rapidity_3122.dat"
# MUSIC input of an event (music_input_mode_2 in Vitor's folders): the first file matching this at the event's top, or one level below
MUSIC_INPUT_PATTERN="*music_input*"
# Keys of the MUSIC inputs that are summarized: the initial condition, then the bullet block (README, "MUSIC bullet settings")
MUSIC_KEYS=(Initial_Distribution_input_filename total_current_energy total_current_momentum_x total_current_momentum_y total_current_momentum_eta_s bullet_size_x bullet_size_y bullet_size_eta_s bullet_position_x bullet_position_y bullet_position_eta_s bullet_initial_time bullet_final_time)
MAX_SAMPLE_DEPTH=3  # samples are looked for down to this depth below a dataset folder
MAX_NAMES=14        # longest list of names printed before it is cut
MAX_RANGES=8        # longest list of event ranges printed before it is cut
MAX_EXAMPLES=5      # events quoted when a check fails for some of them
MAX_LINE_CHARS=150  # lines quoted from parameter files are cut at this length

# Numbers are parsed in the C locale: under a locale with a decimal comma (pt_BR, say), awk reads 4.5e-01 as 4 and every range is wrong
export LC_ALL=C

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

# Prints the header block from "# Usage:" down to the closing line of hashes (same as runHydroTruthLab.sh)
print_help() {
  awk '/^# Usage:/{show=1} show && /^##########/{exit} show' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

DATASETS=()
PAIRS=() # SAMPLE NOBULLET SAMPLE NOBULLET ...
while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      print_help
      ;;
    --base|-b|--iss|-i)
      if [[ $# -lt 2 ]]; then
        echo "Option $1 needs a value. Use --help for usage." >&2
        exit 1
      fi
      if [[ "$1" == "--base" || "$1" == "-b" ]]; then BASE="$2"; else ISS_NAME="$2"; fi
      shift 2
      ;;
    --pair|-p)
      if [[ $# -lt 3 ]]; then
        echo "Option $1 needs two samples. Use --help for usage." >&2
        exit 1
      fi
      PAIRS+=("$2" "$3")
      shift 3
      ;;
    -*)
      echo "Unknown option: $1. Use --help for usage." >&2
      exit 1
      ;;
    *)
      DATASETS+=("$1")
      shift
      ;;
  esac
done

# ==============================================================================
# HELPERS
# ==============================================================================

shopt -s nullglob # a glob that matches nothing expands to nothing, not to itself
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/surveyHydroDatasets.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

# Prints its arguments separated by spaces, cut after MAX_NAMES
join_capped() {
  local count=$#
  local shown=("${@:1:MAX_NAMES}")
  local text="${shown[*]}"
  if (( count > MAX_NAMES )); then
    text+=" (+$(( count - MAX_NAMES )) more)"
  fi
  echo "$text"
}

# Reads integers, one per line, and prints them as ranges ("0-99, 101-249"), cut after MAX_RANGES ranges
compress_ids() {
  sort -n | awk -v maxRanges="$MAX_RANGES" '
    NR == 1 { first = $1; last = $1; next }
    $1 == last + 1 { last = $1; next }
    { ranges[++n] = (first == last ? first : first "-" last); first = $1; last = $1 }
    END {
      if (NR == 0) { print "none"; exit }
      ranges[++n] = (first == last ? first : first "-" last)
      text = ""
      for (i = 1; i <= n && i <= maxRanges; ++i)
        text = text (i > 1 ? ", " : "") ranges[i]
      if (n > maxRanges)
        text = text ", (+" n - maxRanges " more ranges)"
      print text
    }'
}

# First two lines of every parameter file (.dat) at the top of a folder, but iSS_parameters.dat (an iSS configuration, not a parameter
# table). A file identical to the one of the same name in the folder given as $3 is only named
describe_parameter_files() {
  local dir="$1" indent="$2" referenceDir="${3:-}" file name firstLine secondLine fields
  for file in "$dir"/*.dat; do
    name=$(basename "$file")
    [ "$name" = "iSS_parameters.dat" ] && continue
    if [ -n "$referenceDir" ] && [ -f "$referenceDir/$name" ] && cmp -s "$file" "$referenceDir/$name"; then
      echo "${indent}${name}: identical to the one at the top of the dataset"
      continue
    fi
    firstLine=$(sed -n 1p "$file" | tr -d '\r' | cut -c1-"$MAX_LINE_CHARS")
    secondLine=$(sed -n 2p "$file" | tr -d '\r' | cut -c1-"$MAX_LINE_CHARS")
    # Fields per data row: the producer's jet-file parser expects five (event, n_random, phi, momentum_x, momentum_y)
    fields=$(tr -d '\r' < "$file" | awk '/^[ \t]*#/ || NF == 0 { next } { ++rows; if (rows == 1 || NF < low) low = NF; if (NF > high) high = NF }
      END { if (rows == 0) print "no data rows"; else if (low == high) printf "%d data rows of %d fields", rows, low; else printf "%d data rows of %d to %d fields", rows, low, high }')
    echo "${indent}${name}: ${fields} | ${firstLine} | ${secondLine}"
  done
}

# Prints the sample folders at or below a folder: one with numbered subfolders (events), or a single-event folder (a MUSIC input of its
# own). Numbered folders are never descended into, and neither is a sample
find_samples() {
  local dir="$1" depth="$2" sub
  for sub in "$dir"/*/; do
    if [[ "$(basename "$sub")" =~ ^[0-9]+$ ]]; then
      echo "$dir"
      return
    fi
  done
  if [ -n "$(find "$dir" -maxdepth 1 -type f -name "$MUSIC_INPUT_PATTERN" -print -quit 2>/dev/null)" ]; then
    echo "$dir"
    return
  fi
  (( depth >= MAX_SAMPLE_DEPTH )) && return
  for sub in "$dir"/*/; do
    [ -L "${sub%/}" ] && continue # links could loop
    find_samples "${sub%/}" $(( depth + 1 ))
  done
}

# Rows, columns and the nodes of the first three columns (y, pT, phi) of an iSS file
grid_summary() {
  awk '
    /^#/ { next }
    NF >= 3 {
      ++rows
      columns = NF
      if (!($1 in yNodes)) { yNodes[$1]; ++nY }
      if (!($2 in ptNodes)) { ptNodes[$2]; ++nPt }
      if (!($3 in phiNodes)) { phiNodes[$3]; ++nPhi }
      if (rows == 1) { yLow = yHigh = $1 + 0; ptLow = ptHigh = $2 + 0; phiLow = phiHigh = $3 + 0 }
      if ($1 + 0 < yLow) yLow = $1 + 0; if ($1 + 0 > yHigh) yHigh = $1 + 0
      if ($2 + 0 < ptLow) ptLow = $2 + 0; if ($2 + 0 > ptHigh) ptHigh = $2 + 0
      if ($3 + 0 < phiLow) phiLow = $3 + 0; if ($3 + 0 > phiHigh) phiHigh = $3 + 0
    }
    END {
      if (rows == 0) { print "no data rows"; exit }
      printf "%d rows, %d columns; y: %d nodes in [%g, %g], pT: %d in [%g, %g], phi: %d in [%g, %g]\n", rows, columns, nY, yLow, yHigh, nPt, ptLow, ptHigh, nPhi, phiLow, phiHigh
    }' "$1"
}

# Summary of the MUSIC keys over the events of a sample, from "event key value" lines ($1) and the number of inputs read ($2): one line per
# key that varies, then the constant keys and the absent ones on one line each
summarize_music() {
  awk -v keys="${MUSIC_KEYS[*]}" -v nInputs="$2" '
    {
      key = $2; value = $3; ++count[key]
      if (!((key SUBSEP value) in seen)) { seen[key, value]; ++distinct[key]; if (!(key in example)) example[key] = value }
      if (value ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/) {
        number = value + 0
        if (!(key in low) || number < low[key]) low[key] = number
        if (!(key in high) || number > high[key]) high[key] = number
      } else {
        isText[key] = 1
      }
      # An initial condition in a folder named after the event (event i reads .../i/...) is the usual pattern: count the exceptions
      if (key == "Initial_Distribution_input_filename" && $1 != "single" && index(value, "/" $1 "/") > 0) ++ownFolder
    }
    END {
      n = split(keys, list, " ")
      constant = ""; absent = ""
      for (i = 1; i <= n; ++i) {
        key = list[i]
        if (!(key in count)) { absent = absent " " key; continue }
        partial = count[key] < nInputs ? sprintf(" (in %d/%d inputs)", count[key], nInputs) : ""
        if (distinct[key] == 1) { constant = constant " " key "=" example[key] partial; continue }
        line = sprintf("    varies: %s, %d distinct", key, distinct[key])
        line = line (key in isText ? ", e.g. " example[key] : sprintf(" in [%g, %g]", low[key], high[key]))
        if (key == "Initial_Distribution_input_filename")
          line = line sprintf("; in a folder named after the event for %d/%d", ownFolder, count[key])
        print line partial
      }
      if (constant != "") print "    constant:" constant
      if (absent != "") print "    absent:" absent
    }' "$1"
}

# Pairing of a no_bullet-like sample ($1) with another sample ($2), from their "event key value" files: the initial condition, event by event,
# compared through the last three components of its path (<folder>/<event>/<file>). Pairs that differ only in the prefix of the path (a folder
# moved between machines) are counted as the same and reported
compare_initial_conditions() {
  awk -v maxExamples="$MAX_EXAMPLES" '
    function tail(path,    parts, n) { n = split(path, parts, "/"); return n >= 3 ? parts[n - 2] "/" parts[n - 1] "/" parts[n] : path }
    FNR == NR { if ($2 == "Initial_Distribution_input_filename") condition[$1] = $3; next }
    $2 == "Initial_Distribution_input_filename" && ($1 in condition) {
      ++shared
      if (tail(condition[$1]) != tail($3)) {
        if (++different <= maxExamples && $1 != "single") examples = examples " " $1
        next
      }
      ++same
      if (condition[$1] != $3 && ++moved == 1) {
        movedExample = substr($3, 1, length($3) - length(tail($3)) - 1) " vs " substr(condition[$1], 1, length(condition[$1]) - length(tail(condition[$1])) - 1)
      }
    }
    END {
      if (shared == 0) { print "no event with an initial condition in both"; exit }
      printf "%d event(s) in both, same initial condition in %d", shared, same
      if (moved > 0) printf " (%d of them under another prefix: %s)", moved, movedExample
      if (different > 0) {
        printf ", different in %d", different
        if (examples != "") printf " (events%s%s)", examples, (different > maxExamples ? " ..." : "") # parentheses: a bare > is a redirect
      }
      printf "\n"
    }' "$1" "$2"
}

# Surveys one sample ($1, labelled $2) of the dataset at $6: events, parameter files, iSS file, grid and MUSIC inputs. Writes the
# "event key value" lines of its MUSIC inputs to $3 and its grid to $4, and reads the grids printed so far from $5 ("label|grid" lines) to
# print a repeated grid as "same as"
survey_sample() {
  local sample="$1" label="$2" musicValues="$3" gridLine="$4" gridsSoFar="$5" root="$6"
  local ids=() entry single=0
  for entry in "$sample"/*/; do
    entry=$(basename "$entry")
    [[ "$entry" =~ ^[0-9]+$ ]] && ids+=("$entry")
  done
  if (( ${#ids[@]} == 0 )); then
    single=1
    ids=("single")
    echo "  sample ${label}: single event"
  else
    mapfile -t ids < <(printf '%s\n' "${ids[@]}" | sort -n)
    echo "  sample ${label}: ${#ids[@]} events, $(printf '%s\n' "${ids[@]}" | compress_ids)"
  fi
  local firstDir="$sample"
  (( single )) || firstDir="$sample/${ids[0]}"
  [ "$sample" != "$root" ] && describe_parameter_files "$sample" "    " "$root" # the top of the dataset is described already

  # iSS file: at results/ of every event, or wherever the first event has it
  local issPath="" found
  if [ -f "$firstDir/results/$ISS_NAME" ]; then
    issPath="results/$ISS_NAME"
  else
    found=$(find "$firstDir" -maxdepth 3 -type f -name "$ISS_NAME" -print -quit 2>/dev/null)
    [ -n "$found" ] && issPath="${found#"$firstDir"/}"
  fi
  if [ -z "$issPath" ]; then
    local resultFiles=("$firstDir"/results/*)
    if (( ${#resultFiles[@]} > 0 )); then
      echo "    iSS file: not in the first event, whose results/ holds $(join_capped "${resultFiles[@]##*/}")"
    else
      echo "    iSS file: not in the first event, which has no results/"
    fi
  else
    local missing=() id dir
    for id in "${ids[@]}"; do
      dir="$sample/$id"
      (( single )) && dir="$sample"
      [ -f "$dir/$issPath" ] || missing+=("$id")
    done
    local otherSmu=() smuFile
    for smuFile in "$(dirname "$firstDir/$issPath")"/Smu*; do
      [ "$(basename "$smuFile")" != "$ISS_NAME" ] && otherSmu+=("$smuFile")
    done
    local line="    iSS file: ${issPath} in $(( ${#ids[@]} - ${#missing[@]} ))/${#ids[@]} events"
    (( ${#missing[@]} > 0 )) && ! (( single )) && line+=", missing in $(printf '%s\n' "${missing[@]}" | compress_ids)"
    (( ${#otherSmu[@]} > 0 )) && line+="; other Smu files next to it: $(join_capped "${otherSmu[@]##*/}")"
    echo "$line"
    local grid previous
    grid=$(grid_summary "$firstDir/$issPath")
    previous=$(awk -F'|' -v grid="$grid" '$2 == grid { print $1; exit }' "$gridsSoFar")
    if [ -n "$previous" ]; then
      echo "    grid: same as ${previous}"
    else
      echo "    grid: ${grid}"
    fi
    echo "${label}|${grid}" > "$gridLine"
  fi

  # MUSIC inputs: the file the first event has, at its top or one level below, read in every event
  local musicPath="" matches=()
  mapfile -t matches < <(find "$firstDir" -maxdepth 1 -type f -name "$MUSIC_INPUT_PATTERN" 2>/dev/null | sort)
  (( ${#matches[@]} == 0 )) && mapfile -t matches < <(find "$firstDir" -mindepth 2 -maxdepth 2 -type f -name "$MUSIC_INPUT_PATTERN" 2>/dev/null | sort)
  if (( ${#matches[@]} == 0 )); then
    echo "    MUSIC input: none in the first event"
    return
  fi
  musicPath="${matches[0]#"$firstDir"/}"
  local nInputs=0 file missingInputs=()
  : > "$musicValues"
  for id in "${ids[@]}"; do
    file="$sample/$id/$musicPath"
    (( single )) && file="$sample/$musicPath"
    if [ ! -f "$file" ]; then
      missingInputs+=("$id")
      continue
    fi
    nInputs=$(( nInputs + 1 ))
    awk -v event="$id" -v keys="${MUSIC_KEYS[*]}" '
      BEGIN { n = split(keys, list, " "); for (i = 1; i <= n; ++i) wanted[list[i]] = 1 }
      ($1 in wanted) && !($1 in done) { done[$1] = 1; value = $2; sub(/#.*/, "", value); print event, $1, value }' "$file" >> "$musicValues"
  done
  local line="    MUSIC input: ${musicPath} in ${nInputs}/${#ids[@]} events"
  if (( ${#matches[@]} > 1 )); then
    local otherInputs=("${matches[@]:1}")
    line+=" (the first event also has $(join_capped "${otherInputs[@]#"$firstDir"/}"))"
  fi
  (( ${#missingInputs[@]} > 0 )) && ! (( single )) && line+=", missing in $(printf '%s\n' "${missingInputs[@]}" | compress_ids)"
  echo "$line"
  summarize_music "$musicValues" "$nInputs"
}

# Surveys one dataset folder ($1, labelled $2): its top, every sample, and the pairings
survey_dataset() {
  local root="$1" name="$2"
  echo
  echo "=== ${name} ==="
  if [ ! -d "$root" ]; then
    echo "  not a folder"
    return
  fi
  local entries=() entry base
  for entry in "$root"/*; do
    base=$(basename "$entry")
    [[ -d "$entry" && "$base" =~ ^[0-9]+$ ]] && continue # events: reported with their sample below
    [ -d "$entry" ] && base+="/"
    entries+=("$base")
  done
  echo "  top: $(join_capped "${entries[@]}")"
  describe_parameter_files "$root" "  "

  local samples=()
  mapfile -t samples < <(find_samples "$root" 0)
  if (( ${#samples[@]} == 0 )); then
    echo "  no samples (no numbered event folders and no MUSIC input down to depth ${MAX_SAMPLE_DEPTH})"
    return
  fi
  local datasetTmp="$TMP_DIR/${name//\//_}"
  mkdir -p "$datasetTmp"
  : > "$datasetTmp/grids"
  local labels=() index label
  for index in "${!samples[@]}"; do
    label="${samples[$index]#"$root"}"
    label="${label#/}"
    [ -z "$label" ] && label="."
    labels+=("$label")
    : > "$datasetTmp/grid$index"
    survey_sample "${samples[$index]}" "$label" "$datasetTmp/music$index" "$datasetTmp/grid$index" "$datasetTmp/grids" "$root"
    cat "$datasetTmp/grid$index" >> "$datasetTmp/grids"
  done

  # Pairing: every no_bullet-like sample against every other sample with the same parent folder
  local other
  for index in "${!samples[@]}"; do
    [[ "$(basename "${samples[$index]}")" == *no_bullet* ]] || continue
    [ -s "$datasetTmp/music$index" ] || continue
    for other in "${!samples[@]}"; do
      [[ "$(basename "${samples[$other]}")" == *no_bullet* ]] && continue
      [ "$(dirname "${samples[$other]}")" = "$(dirname "${samples[$index]}")" ] || continue
      [ -s "$datasetTmp/music$other" ] || continue
      echo "  pairing ${labels[$other]} vs ${labels[$index]}: $(compare_initial_conditions "$datasetTmp/music$index" "$datasetTmp/music$other")"
    done
  done
}

# ==============================================================================
# MAIN
# ==============================================================================

if [ ! -d "$BASE" ]; then
  echo "Error: ${BASE} is not a folder" >&2
  exit 1
fi
# Everything runs inside the temporary folder, so that not even a stray relative write can land among the data or in the caller's folder
BASE=$(cd "$BASE" && pwd)
cd "$TMP_DIR" || exit 1
if (( ${#DATASETS[@]} == 0 && ${#PAIRS[@]} == 0 )); then
  for entry in "$BASE"/*/; do
    DATASETS+=("$(basename "$entry")")
  done
fi
echo "Survey of ${BASE}, $(date '+%Y-%m-%d %H:%M'), iSS file ${ISS_NAME}"
for dataset in "${DATASETS[@]}"; do
  survey_dataset "$BASE/$dataset" "$dataset"
done
# Explicit pairs, also across dataset folders: both samples surveyed, then compared
for (( pairIndex = 0; pairIndex < ${#PAIRS[@]}; pairIndex += 2 )); do
  pairDir="$TMP_DIR/pair$(( pairIndex / 2 ))"
  mkdir -p "$pairDir"
  : > "$pairDir/grids"
  echo
  echo "=== pair ${PAIRS[$pairIndex]} vs ${PAIRS[$pairIndex + 1]} ==="
  for side in 0 1; do
    samplePath="${PAIRS[$pairIndex + side]}"
    [[ "$samplePath" = /* ]] || samplePath="$BASE/$samplePath"
    if [ ! -d "$samplePath" ]; then
      echo "  ${PAIRS[$pairIndex + side]}: not a folder"
      continue 2
    fi
    : > "$pairDir/grid$side"
    survey_sample "$samplePath" "${PAIRS[$pairIndex + side]}" "$pairDir/music$side" "$pairDir/grid$side" "$pairDir/grids" "$(dirname "$samplePath")"
    cat "$pairDir/grid$side" >> "$pairDir/grids"
  done
  if [ -s "$pairDir/music0" ] && [ -s "$pairDir/music1" ]; then
    echo "  pairing: $(compare_initial_conditions "$pairDir/music1" "$pairDir/music0")"
  else
    echo "  pairing: no MUSIC inputs to compare"
  fi
done