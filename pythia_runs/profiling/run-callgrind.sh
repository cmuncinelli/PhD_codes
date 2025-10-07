#!/usr/bin/env bash
EXE=../PythiaGenMin.exe
ARGS="\"/storage3/cicero/pythia_data\" \"/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg\" 1e4 324"
OUT=callgrind.out
# Running single core for valgrind to work safely!
# By the way, should compile with:
# g++ -O2 -std=c++17 -march=native -flto=2 -pipe -fopenmp -DCGAL_NO_INTERVAL_TEST PythiaGenMin.cxx -o PythiaGenMin.exe -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8 -I$FASTJET/include -L$FASTJET/lib -lfastjet `root-config --cflags --libs`

valgrind --tool=callgrind --dump-instr=yes --collect-jumps=yes \
  --callgrind-out-file=$OUT $EXE $ARGS

echo "Callgrind output in $OUT. Visualize with: kcachegrind $OUT"