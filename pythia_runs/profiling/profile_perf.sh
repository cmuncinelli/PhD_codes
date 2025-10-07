#!/usr/bin/env bash
# # This is meant to run outside of the profiling/ folder, in the pythia_runs folder, so that PythiaGenMin.exe is in the current folder!
# EXE=./PythiaGenMin.exe
# ARGS="\"/storage3/cicero/pythia_data\" \"/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg\" 1e4 324"
# OUT=perf.data
# FOLDED=perf-folded.txt
# SVG=perf-report.svg

# echo "Running perf record (this will run your program once)..."
# # Record stack samples, capture user-space callstacks
# perf record -g -- $EXE $ARGS -o /dev/null

# echo "Generating report..."
# # interactive: perf report
# # folded stack format for flamegraph tools:
# perf script > perf.unfolded
# # Create folded stacks (needs stackcollapse-perf.pl from FlameGraph repo)
# # If you have FlameGraph tools, do:
# if command -v stackcollapse-perf.pl >/dev/null 2>&1; then
#   stackcollapse-perf.pl perf.unfolded > $FOLDED
#   echo "Folded stacks written to $FOLDED (use FlameGraph to render)."
# else
#   echo "stackcollapse-perf.pl not found — you'll need it (FlameGraph repository) to make flamegraphs."
# fi

# echo "You can inspect with: perf report --stdio"


# Actually a simpler version of this:
perf record -g \
./PythiaGenMin.exe "/storage3/cicero/pythia_data" "/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg" 1e5 324

echo "Now run \"perf report\" to see the output!"
echo "Or even something like:"
echo "perf report -i perf.data --stdio --sort comm,dso,symbol | head -n 40"
echo "perf report -i perf.data --stdio --sort overhead | head -n 20"
echo "perf report -i perf.data --stdio > report.txt"
echo "for shorter output, or output that is sorted by CPU/usage intensive top 40 or top 20 functions"
echo "(Or just run the extract_top_functions_from_profile_perf.sh script!)"