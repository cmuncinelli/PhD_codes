#!/usr/bin/env bash
perf stat -e cache-references,cache-misses -o ./profiling/perf-summary1.txt \
./PythiaGenMin.exe "/storage3/cicero/pythia_data" "/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg" 1e5 324

# Another version, with some useful flags:
perf stat -e \
instructions,cycles,cache-references,cache-misses,branch-instructions,branch-misses,context-switches,cpu-migrations,page-faults,minor-faults,major-faults,L1-dcache-loads,L1-dcache-load-misses,L1-icache-loads \
-o ./profiling/perf-summary2.txt \
./PythiaGenMin.exe "/storage3/cicero/pythia_data" "/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg" 1e5 324

# Combining the two results into a single .txt file:
cat ./profiling/perf-summary1.txt ./profiling/perf-summary2.txt > ./profiling/perf-summary-merged.txt

# Deleting the unmerged files:
rm ./profiling/perf-summary1.txt ./profiling/perf-summary2.txt