#!/bin/bash
run_name="Asymmetry"
LOGFILE="log-${run_name}.txt"
DIR_THIS=$PWD

FileIn="$1"
JSON="$2"
OPTION="-b --configuration json://${JSON}"

# New tests using tasks that should be able to run locally the TOF data:
o2-analysis-lf-strangenesstofpid ${OPTION} |
o2-analysis-lf-v0mlscoresconverter ${OPTION} |
o2-analysis-lf-asymmetric-rapidity-test ${OPTION} --aod-file $FileIn 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
    echo "No problems!"
    mkdir -p "${DIR_THIS}/results/"
    mv AnalysisResults.root "${DIR_THIS}/results/AnalysisResults.root"
    mv dpl-config.json "${DIR_THIS}/results/${run_name}.json"
else
    echo "Error: Exit code ${rc}"
    echo "Check the log file ${LOGFILE}"
    exit ${rc}
fi

# Will run like this:
# ./run_test.sh "/home/cicero/Datasets/asymmetric_rapidity_tests/AO2D_1.root" "dpl-config-Asymmetry_test.json"
    # Or for a list of files inside a .txt:
# ./run_test.sh "@input_data.txt" "dpl-config-Asymmetry_test.json"

    # Generating dpl-config.json:
# o2-analysis-lf-lambdainvmasstest | o2-analysis-lf-strangenesstofpid | o2-analysis-lf-v0mlscoresconverter --aod-file /home/cicero/Datasets/AO2D_strangeness_data_sets/AO2D_LHC23_pass4_Thin.root -b

# Code runs as:
# o2-analysis-lf-asymmetric-rapidity-test | o2-analysis-lf-strangenesstofpid | o2-analysis-lf-v0mlscoresconverter --aod-file AO2D_LHC23_pass4_Thin.root -b