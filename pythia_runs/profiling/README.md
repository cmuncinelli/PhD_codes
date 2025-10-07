The [**profile_perf.sh**](profile_perf.sh) and [**run-callgrind.sh**](run-callgrind.sh) codes are meant for profiling the code and getting the regions where it spends most of the time running.

You should probably just use [**simpler_perf.sh**](simpler_perf.sh) for simpler checks and a general feel of what the code is doing.

### Usage for these profiling codes is:
(chmod +x profile_perf.sh)
./profile_perf.sh
perf report --stdio    # open interactive textual report

(install valgrind kcachegrind via apt first)
(chmod +x run-callgrind.sh)
./run-callgrind.sh

Open the callgrind.out with kcachegrind or callgrind_annotate for textual reports.

-- Or a simpler command: perf stat -- ./PythiaGenMin.exe "/storage3/cicero/pythia_data" "/home/users/cicerodm/PhD_codes/pythia_runs/input_cards_ALICE/pythia8_pp_136tev_AllowLambdaDecays.cfg" 1e4 324
-- To run it, check the paranoid level with sysctl -n kernel.perf_event_paranoid (SHOULD BE 4 !), then do sudo sysctl -w kernel.perf_event_paranoid=0, run your code in your regular user, and then reset with sudo sysctl -w kernel.perf_event_paranoid=4

Notes:
1. perf record -F 99 -g samples at 99 Hz -- adjust (-F) if you want denser sampling.
2. If perf cannot record due to permissions, run with sudo